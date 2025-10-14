# Genome Analyzer Pipeline

## 1. 소개 (Introduction)

Genome Analyzer Pipeline은 박테리아 유전체 서열의 분자 역학(Molecular Epidemiology) 분석을 자동화하기 위해 설계된 Python 기반 커맨드 라인 도구입니다. FASTA 파일을 입력받아 MLST, 항생제 내성 유전자 탐지 등 일련의 생물정보학 분석을 수행하고, 그 결과를 하나의 포괄적인 요약 보고서로 통합합니다.

이 파이프라인의 주요 목표는 반복적인 분석 작업을 자동화하여 연구 효율성을 높이고, 박테리아 유전체에 대한 빠르고 재현 가능한(reproducible) 인사이트를 제공하는 것입니다.

---

## 2. 주요 기능 (Key Features)

- **자동 MLST 분석**: 7-gene MLST 스킴을 기반으로 housekeeping gene을 자동으로 식별하고, allele 번호를 결정하며, 최종 Sequence Type (ST)을 부여합니다.
- **동적 종(Species) 식별**: 입력된 유전체의 FASTA 헤더에서 박테리아 종을 지능적으로 감지하여, 해당 종에 적합한 MLST 데이터베이스를 동적으로 선택하여 분석합니다.
- **포괄적인 유전자 탐지**: 다음을 포함한 광범위한 중요 유전자를 식별합니다:
    - 항생제 내성 (AMR, Antimicrobial Resistance) 유전자
    - 플라스미드 복제원 (Plasmid Replicons)
    - 이동성 유전 인자 (MGEs, Mobile Genetic Elements)
- **고성능 비동기(Asynchronous) 처리**: Python의 `asyncio` 라이브러리를 활용하여 여러 BLAST 분석을 동시에 실행함으로써 총 분석 시간을 크게 단축합니다.
- **모듈화 및 확장 가능한 아키텍처**: 코드가 기능별 모듈(module)로 명확히 구성되어 있어, 새로운 분석 유형을 추가하거나 기존 로직을 수정하기 용이합니다.
- **상세 로깅**: 전용 로거(logger)가 각 분석 모듈의 단계별 출력을 캡처하여 디버깅을 용이하게 하고 투명한 감사 추적(audit trail)을 제공합니다. `--verbose` 플래그는 실시간 콘솔 피드백도 제공합니다.
- **견고한 오류 처리**: 외부 의존성(NCBI BLAST+)을 확인하고, 사용자 안내를 위해 명확하고 유용한 오류 메시지를 제공합니다.

---

## 3. 프로젝트 구조 (Project Structure)

프로젝트는 소스 코드, 데이터베이스, 입출력 데이터를 위한 별도의 디렉토리로 구성되어 명확성과 유지보수성을 높입니다.

```
genome_analyzer/
├── src/                    # 1. 메인 소스 코드
│   ├── main.py             #    - CLI 진입점 및 argument 파싱
│   ├── config.py           #    - 경로 및 분석을 위한 중앙 설정
│   ├── logger.py           #    - 단계별 디버그 로거
│   ├── analysis/           #    - 핵심 분석 로직 모듈
│   │   ├── manager.py      #    - 전체 분석 파이프라인 조율
│   │   ├── blast_runner.py #    - BLAST 커맨드를 비동기적으로 실행
│   │   └── utils.py        #    - 헬퍼 함수 (e.g., 의존성 체크, 종 식별)
│   └── reporting/          #    - 리포트 생성 모듈
│       └── reporter.py     #    - 결과를 최종 텍스트 리포트로 컴파일
│
├── database/               # 2. 모든 분석 데이터베이스의 루트 디렉토리
│   ├── MLST_DB/            #    - MLST 스킴 (종별로 하위 폴더 구성)
│   │   └── klebsiella/
│   ├── resfinder_db/       #    - 항생제 내성 유전자
│   ├── plasmidfinder_db/   #    - 플라스미드 복제원 서열
│   └── mefinder_db/        #    - 이동성 유전 인자
│
├── genome/                 # 3. 입력 유전체 파일 위치
│   └── GCF_000523395.1.fna
│
├── analysis_results/       # 4. 리포트 및 결과물의 기본 출력 디렉토리
│
└── logs/                   # 5. 디버깅 및 감사 추적을 위한 상세 로그
```

---

## 4. 코드 로직 및 모듈 설명

파이프라인의 로직은 모듈식이며, 각 컴포넌트가 워크플로우의 특정 부분을 담당합니다.

- **`main.py` (Entry Point)**
  - `argparse`를 사용하여 커맨드 라인 인터페이스(`genome_file`, `--output`, `--verbose`)를 정의합니다.
  - 제공된 인자로 `AnalysisManager`를 인스턴스화합니다.
  - `asyncio.run()`을 사용하여 비동기 파이프라인을 시작합니다.

- **`config.py` (Central Configuration)**
  - `DATABASE_ROOT`, `DEFAULT_RESULTS_DIR`와 같은 핵심 경로를 정의합니다.
  - 파이프라인의 **컨트롤 센터** 역할을 하는 `ANALYSES_TO_RUN` 딕셔너리를 포함합니다. 이 딕셔너리는 데이터베이스 폴더와 분석 이름을 매핑하여, 파이프라인이 데이터 기반으로 동작하고 쉽게 확장될 수 있도록 합니다.

- **`analysis/manager.py` (The Orchestrator)**
  - `AnalysisManager` 클래스는 전체 워크플로우를 총괄합니다. `run_pipeline` 메소드는 다음 순서로 작업을 실행합니다:
    1.  **Setup**: 출력 디렉토리를 생성하고 `utils.check_dependencies()`를 통해 사전 검사를 수행합니다.
    2.  **Genome DB Creation**: 효율적인 검색을 위해 입력 유전체 파일로부터 로컬 BLAST 데이터베이스를 생성합니다.
    3.  **Species Identification**: `utils.setup_mlst_parameters()`를 호출하여 종(species)을 결정하고 MLST 관련 데이터를 준비합니다.
    4.  **Concurrent Analysis**: 핵심 성능 기능입니다. MLST 워크플로우와 `config.py`에 정의된 다른 모든 분석을 위한 `asyncio` task 목록을 생성합니다. `asyncio.gather()`는 BLAST 위주의 모든 task를 동시에 실행합니다.
    5.  **Reporting**: 모든 분석이 완료되면 `reporter.create_final_report()`를 호출하여 요약 리포트를 생성합니다.
    6.  **Cleanup**: 분석 중 사용된 임시 디렉토리를 제거합니다.

- **`analysis/blast_runner.py` (The Asynchronous Worker)**
  - NCBI BLAST+ 커맨드 라인 툴(`makeblastdb`, `blastn`)에 대한 비동기 래퍼(wrapper)를 제공합니다.
  - `asyncio.create_subprocess_exec`를 사용하여 셸 커맨드를 non-blocking 방식으로 실행하고, `stdout`과 `stderr`를 캡처하여 견고한 오류 처리를 구현합니다. 이를 통해 여러 개의 긴 BLAST 작업을 병렬로 실행할 수 있습니다.

- **`analysis/utils.py` (Helper Functions)**
  - `check_dependencies()`: 필수 BLAST 툴들이 시스템의 PATH에 설치되어 있는지 확인합니다.
  - `setup_mlst_parameters()`: 입력 FASTA 헤더를 읽어 종(e.g., "klebsiella")을 자동으로 식별하는 핵심 함수입니다. 이후 해당 MLST 데이터베이스 디렉토리와 프로파일 파일을 동적으로 찾아 파이프라인이 다른 유기체에도 적용될 수 있도록 합니다.

- **`reporting/reporter.py` (The Scribe)**
  - `AnalysisManager`가 수집한 구조화된 Python 딕셔너리 형태의 결과를 입력받습니다.
  - 각 분석의 BLAST 결과 파일(`.tsv`)을 파싱합니다.
  - ST, allele profile, AMR 유전자, plasmid replicon 등 모든 정보를 사람이 읽기 쉬운 `Final_ME_Report.txt` 파일로 포맷팅합니다.

- **`logger.py` (The Archivist)**
  - `logs/` 디렉토리에 상세한 단계별 로그 파일을 기록하는 `Logger` 클래스를 제공합니다.
  - 각 로그 파일은 `2025-10-14_MLST_4_Extracted_Genes_Content_1.fasta`와 같이 체계적인 이름으로 저장되어, 중간 데이터와 커맨드 출력을 명확하게 추적할 수 있게 해주므로 디버깅에 매우 유용합니다.

---

## 5. 시스템 요구사항 (System Requirements)

- **Python**: 3.9 이상
- **NCBI BLAST+ Suite**: `blastn`, `makeblastdb`, `blastdbcmd`가 설치되어 시스템의 PATH에서 접근 가능해야 합니다.
- **Python Libraries**:
  - `pandas`
  - `biopython`

---

## 6. 설치 및 설정 (Installation and Setup)

**1. Repository 클론**
```bash
git clone <repository_url>
cd genome_analyzer
```

**2. Python 라이브러리 설치**
```bash
pip install pandas biopython
```

**3. 데이터베이스 설정**
프로젝트 루트에 `database` 디렉토리를 생성합니다. 내부에 다음 구조에 따라 필요한 데이터베이스 파일을 배치합니다. MLST 데이터베이스는 **반드시** 종의 이름을 딴 하위 디렉토리(e.g., `klebsiella`)를 가져야 합니다.

```
database/
├── MLST_DB/
│   └── klebsiella/         # 종(species)별 폴더
│       ├── gapA.tfa
│       ├── infB.tfa
│       └── klebsiella.txt  # 프로파일 정의 파일
│
├── resfinder_db/           # ResFinder 데이터베이스 파일 (*.fsa)
│   └── all.fsa
│
├── plasmidfinder_db/       # PlasmidFinder 데이터베이스 파일 (*.fsa)
│   └── enterobacteriales.fsa
│
└── mefinder_db/            # MGE 데이터베이스 파일 (*.fasta)
    └── MGEdb_cds.fasta
```

---

## 7. 사용 예제 (Usage Examples)

파이프라인은 `src/main.py`를 통해 실행됩니다.

**1. 도움말 정보 보기**
```bash
python src/main.py --help
```

**2. 기본 실행**
유전체 파일에 대한 분석을 실행합니다. 결과는 기본 `analysis_results/` 디렉토리에 저장됩니다.

```bash
python src/main.py path/to/your/genome.fna
```

**3. 고급 실행**
사용자 지정 출력 디렉토리(`-o`)를 지정하고, 콘솔에 상세 로깅을 활성화(`-v`)합니다.

```bash
python src/main.py path/to/your/genome.fna -o custom_results -v
```

---

## 8. 출력 설명 (Output Description)

완료 시 파이프라인은 지정된 출력 디렉토리에 다음을 생성합니다:

- **`Final_ME_Report.txt`**: MLST, AMR 유전자, MGE를 포함한 모든 분석 결과의 통합된, 사람이 읽기 쉬운 요약 리포트입니다.
- **분석별 폴더**: 각 분석을 위한 디렉토리(e.g., `Antimicrobial_Resistance/`, `Plasmid_Replicons/`)로, 원시 BLAST 결과(`blast_results.tsv`)와 검색에 사용된 쿼리 파일을 포함합니다.
- **`logs/` 디렉토리**: 프로젝트 루트에 위치하며, 각 단계의 상세 로그를 포함하여 문제 해결이나 분석 과정 감사에 유용합니다.

---

## 9. 확장 가이드 (Expansion Guide)

모듈식 설계 덕분에 새로운 BLAST 기반 분석을 추가하는 것은 매우 간단합니다.

**1단계: 데이터베이스 추가**
새로운 FASTA 데이터베이스 파일을 `database/` 디렉토리 안의 전용 폴더에 배치합니다.
```
database/
└── virulence_db/           # 새로운 데이터베이스 폴더
    └── vfdb_core.fasta
```

**2단계: 설정 업데이트**
`src/config.py`를 열고 `ANALYSES_TO_RUN` 딕셔너리에 새 항목을 추가합니다. `key`는 데이터베이스 폴더 이름이고, `value`는 출력 폴더 및 리포트 섹션에 사용할 이름입니다.

```python
# in src/config.py
ANALYSES_TO_RUN = {
    "MLST_DB": "MLST",
    "resfinder_db": "Antimicrobial_Resistance",
    "plasmidfinder_db": "Plasmid_Replicons",
    "mefinder_db": "Mobile_Genetic_Elements",
    "virulence_db": "Virulence_Factors", # <-- 여기에 새 분석을 추가
}
```

이제 모든 준비가 끝났습니다. `AnalysisManager`는 새 항목을 자동으로 감지하여 다른 분석과 동시에 실행하며, `reporter`는 최종 리포트에 해당 결과를 포함시킵니다. 다른 코드를 변경할 필요가 없습니다.

---

## 10. 특별 분석 프로세스 추가하기 (Adding a Special Analysis Process)

9번 가이드에서 설명한 "표준(Standard)" 분석 추가 방법 외에도, 파이프라인은 MLST와 같이 여러 단계의 복잡한 로직을 갖는 "특별(Special)" 분석을 추가할 수 있도록 설계되었습니다. MLST가 바로 그 예시입니다.

현재 아키텍처에서 새로운 특별 분석(예: Serotyping)을 추가하는 방법은 다음과 같습니다.

**1단계: `config.py`에 분석 정의**

먼저 `src/config.py`의 `ANALYSES_TO_RUN` 딕셔너리에 새로운 분석을 추가합니다. 이 과정은 표준 분석을 추가할 때와 동일합니다.

```python
# in src/config.py
ANALYSES_TO_RUN = {
    "MLST_DB": "MLST",
    "serotype_db": "Serotyping", # <-- 새로운 특별 분석 추가
    ...
}
```

**2단계: `AnalysisManager`에 전용 워크플로우 메소드 생성**

`src/analysis/manager.py`의 `AnalysisManager` 클래스 내에, 새로운 특별 분석을 위한 전용 `async` 메소드를 작성합니다. 이 메소드는 `_run_mlst_workflow`를 모델로 삼아, 해당 분석에 필요한 모든 커스텀 로직(예: 여러 BLAST 실행, 결과 파싱, 외부 툴 호출 등)을 포함해야 합니다.

```python
# in src/analysis/manager.py

class AnalysisManager:
    # ... 기존 코드 ...

    async def _run_mlst_workflow(self, genome_db_path: Path, mlst_params: dict):
        # ... MLST 로직 ...

    # --- 새로운 특별 분석 메소드 추가 --- 
    async def _run_serotyping_workflow(self, genome_db_path: Path, serotype_params: dict):
        self._log("Starting Serotyping workflow...")
        # 1. Serotyping에 필요한 첫 BLAST 실행
        # 2. 결과 파싱
        # 3. 파싱된 결과를 바탕으로 두 번째 BLAST 또는 다른 로직 수행
        # 4. 최종 결과를 self.results_data에 저장
        self.results_data['serotype'] = { "Type": "ST1a" } # 예시 결과
        self._log("Serotyping workflow completed.")
```

**3단계: `run_pipeline` 로직에 새 워크플로우 연결**

마지막으로, `AnalysisManager`의 `run_pipeline` 메소드 안에서 `ANALYSES_TO_RUN`을 순회하는 메인 루프를 찾아, 새로운 분석을 호출하는 `elif` 구문을 추가합니다.

```python
# in run_pipeline() method of AnalysisManager

        tasks = []
        for db_folder, analysis_name in ANALYSES_TO_RUN.items():
            if analysis_name == "MLST":
                tasks.append(self._run_mlst_workflow(genome_db_path, mlst_params))

            # --- 여기에 새로운 elif 구문 추가 ---
            elif analysis_name == "Serotyping":
                # serotype_params = utils.setup_serotype_parameters(...) # 필요 시 파라미터 준비
                tasks.append(self._run_serotyping_workflow(genome_db_path, {}))
            # -------------------------------------

            else:
                tasks.append(self._run_other_analysis(db_folder, analysis_name, genome_db_path))
        
        await asyncio.gather(*tasks)
```

이 3단계를 통해, 현재 아키텍처 내에서 단순한 BLAST 검색 이상의 복잡한 로직을 가진 새로운 분석 파이프라인을 효과적으로 추가하고 관리할 수 있습니다.
