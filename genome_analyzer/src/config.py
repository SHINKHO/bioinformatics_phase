# 파이프라인 중앙 설정.
from pathlib import Path

# 프로젝트 루트 디렉토리.
PROJECT_ROOT = Path.cwd()

# 주요 데이터 및 결과 디렉토리.
DATABASE_ROOT = PROJECT_ROOT / "database"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "analysis_results"

# 입력 게놈으로 생성된 BLAST 데이터베이스 저장 디렉토리.
BLAST_DB_DIR = PROJECT_ROOT / "blast_db_output"
BLAST_DB_DIR.mkdir(exist_ok=True)

# 실행할 분석 목록을 정의하는 메인 컨트롤 딕셔너리.
# Key: database/ 내의 DB 폴더 이름
# Value: 출력 폴더 및 보고서에 사용될 분석 이름
ANALYSES_TO_RUN = {
    # 특별 분석 (MLSTHandler가 처리)
    "MLST_DB": "MLST",
    # "Pathogenfinder": "Pathogen_Finder2",
    # 표준 분석 (StandardAnalysisHandler가 처리)
    "resfinder_db": "Antimicrobial_Resistance",
    "plasmidfinder_db": "Plasmid_Replicons",
    "mefinder_db": "Mobile_Genetic_Elements",
}
