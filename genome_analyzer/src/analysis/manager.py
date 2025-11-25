# 메인 분석 오케스트레이터
import asyncio
import time
from pathlib import Path
import shutil

# Local (sibling) module imports
from . import utils
from . import blast_runner
from .handler import (AnalysisContext, MLSTHandler, PathogenFinder2Handler, 
                          StandardAnalysisHandler, AMRHandler)

# Project-level module imports
from reporting import reporter
from config import ANALYSES_TO_RUN, BLAST_DB_DIR
from logger import Logger

class AnalysisManager:
    # 유전체 분석 워크플로우 전체를 관리.
    def __init__(self, genome_file: Path, results_dir: Path, verbose: bool = False):
        # 분석 관리자 초기화. In: genome_file(Path), results_dir(Path), verbose(bool) / Out: None
        self.genome_file = genome_file
        self.base_results_dir = results_dir # Store base results dir
        self.results_dir = results_dir # To be updated
        self.temp_dir = results_dir / "temp" # To be updated
        self.verbose = verbose
        self.results_data = {}
        
        # Set up logging
        self.base_logs_dir = Path.cwd() / "logs"
        self.logs_dir = self.base_logs_dir
        self.logger = Logger(self.logs_dir)

    def _log(self, message: str, level: str = "INFO"):
        # verbose 모드일 때 로그 출력. In: message(str), level(str) / Out: None
        if self.verbose:
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp} - {level}] {message}")

    async def run_pipeline(self):
        # 전체 분석 파이프라인 실행. In: None / Out: None
        start_time = time.time()
        print("===== Genome Analysis Pipeline Start =====")
        try:
            # 1. 사전 확인 및 설정
            self._log("Step 1: Pre-flight checks and setup.")
            self.logger.log_step("Pipeline", "1_Pre-flight_Checks", "Starting pre-flight checks and setup.")
            utils.check_dependencies()

            # 2. 종 식별 및 경로 설정
            self._log("Step 2: Identifying species and setting up paths.")
            mlst_params = utils.setup_mlst_parameters(self.genome_file, self.logger)
            self.results_data['mlst_params'] = mlst_params
            genome_id = mlst_params['genome_id']
            species = mlst_params['species']
            self._log(f"Species '{species}' identified for MLST from folder structure.")
            self.logger.log_step("Pipeline", "2_Species_Identification", f"Species '{species}' identified from folder structure.")

            # genome ID와 species에 따라 출력 디렉토리 정의
            self.results_dir = self.base_results_dir / genome_id / species
            self.temp_dir = self.results_dir / "temp"
            self.logs_dir = self.base_logs_dir / genome_id / species
            blast_db_dir = BLAST_DB_DIR / genome_id / species

            # 디렉토리 생성
            self.results_dir.mkdir(parents=True, exist_ok=True)
            self.temp_dir.mkdir(exist_ok=True)
            blast_db_dir.mkdir(parents=True, exist_ok=True)
            
            # 새 경로로 로거 다시 초기화
            self.logger = Logger(self.logs_dir)

            # 3. 입력 게놈에 대한 BLAST DB 생성
            self._log("Step 3: Creating BLAST database for the input genome.")
            self.logger.log_step("Pipeline", "3_Create_BLAST_DB", "Creating BLAST database for the input genome.")
            genome_db_path = await blast_runner.create_blast_db_async(self.genome_file, blast_db_dir)
            self._log(f"Genome BLAST DB created at '{genome_db_path}'.")
            self.logger.log_step("Pipeline", "4_BLAST_DB_Created", f"Genome BLAST DB created at '{genome_db_path}'.")

            # 4. 모든 분석 작업 동시 실행
            self._log("Step 4: Running all analysis tasks concurrently.")
            self.logger.log_step("Pipeline", "5_Run_Concurrent_Analyses", "Running all analysis tasks concurrently.")
            
            # 핸들러에 전달할 컨텍스트 준비
            context = AnalysisContext(
                genome_db_path=genome_db_path,
                results_dir=self.results_dir,
                temp_dir=self.temp_dir,
                logger=self.logger,
                verbose=self.verbose,
                results_data=self.results_data,
                genome_id=genome_id,
                species=species
            )

            # 핸들러 체인 구성: MLST -> AMR -> PathogenFinder2 -> Standard
            standard_handler = StandardAnalysisHandler(context)
            pathogen_handler = PathogenFinder2Handler(context)
            amr_handler = AMRHandler(context)
            analysis_chain = MLSTHandler(context)
            analysis_chain.set_next(amr_handler).set_next(pathogen_handler).set_next(standard_handler)

            # 핸들러 체인에 모든 분석 전달
            tasks = []
            for db_folder, analysis_name in ANALYSES_TO_RUN.items():
                # 분석별 파라미터 전달
                if analysis_name == "MLST":
                    params = mlst_params
                elif analysis_name == "Pathogen_Finder2":
                    params = {
                        "database_dir": str(Path.cwd() / "database" / "Pathogenfinder"),
                        "output_dir": str(self.results_dir / "Pathogen_Finder2"),
                        "genome_path": str(genome_db_path)
                    }
                else:
                    params = {}
                
                # handle 메소드는 실행 준비된 asyncio.Task를 반환
                task = await analysis_chain.handle(
                    analysis_name=analysis_name,
                    db_folder=db_folder,
                    params=params
                )
                if task:
                    tasks.append(task)
            
            # 생성된 모든 태스크를 동시 실행
            await asyncio.gather(*tasks)
            self._log("All analysis tasks completed.")
            self.logger.log_step("Pipeline", "6_Concurrent_Analyses_Complete", "All analysis tasks completed.")

            # 5. 최종 보고서 생성
            self._log("Step 5: Generating final report.")
            self.logger.log_step("Pipeline", "7_Generate_Report", "Generating final report.")
            genome_name = utils.get_genome_name(self.genome_file)
            reporter.create_final_report(self.results_data, self.results_dir, genome_name)
            self._log("Final report generated.")
            self.logger.log_step("Pipeline", "8_Report_Generated", "Final report generated.")

        except (ValueError, FileNotFoundError, RuntimeError, Exception) as e:
            # 오류 처리
            print(f"\n❌ PIPELINE FAILED: An error occurred.\n  -> {e}")
            self.logger.log_step("Pipeline", "9_Pipeline_Failed", f"PIPELINE FAILED: An error occurred.\n  -> {e}")
        
        finally:
            # 6. 임시 파일 정리
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
        
        end_time = time.time()
        print(f"\n🎉 Analysis complete in {end_time - start_time:.2f} seconds.")
        print(f"   Results are in '{self.results_dir}'")
        print("==========================================")
