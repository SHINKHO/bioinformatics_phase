# Genome Analyzer Pipeline

## 1. 소개 (Introduction)

Genome Analyzer Pipeline은 박테리아 유전체 서열의 분자 역학(Molecular Epidemiology) 분석을 자동화하기 위해 설계된 Python 기반 커맨드 라인 도구입니다. FASTA 파일을 입력받아 MLST, 항생제 내성 유전자 탐지 등 일련의 생물정보학 분석을 수행하고, 그 결과를 하나의 포괄적인 요약 보고서로 통합합니다.

이 파이프라인의 주요 목표는 반복적인 분석 작업을 자동화하여 연구 효율성을 높이고, 박테리아 유전체에 대한 빠르고 재현 가능한(reproducible) 인사이트를 제공하는 것입니다.

---

## 2. 주요 기능 (Key Features)

- **자동 MLST 분석**: 7-gene MLST 스킴을 기반으로 housekeeping gene을 자동으로 식별하고, allele 번호를 결정하며, 최종 Sequence Type (ST)을 부여합니다.
- **PathogenFinder2 통합**: 병원체 식별 및 특성화를 위한 PathogenFinder2 분석을 지원하여 병원체 유전체의 종 분류와 위험도 평가를 수행합니다.
- **동적 종(Species) 식별**: 입력된 유전체의 FASTA 헤더에서 박테리아 종을 지능적으로 감지하여, 해당 종에 적합한 MLST 데이터베이스를 동적으로 선택하여 분석합니다.
- **지능적인 샘플 ID 추출 (Intelligent Sample ID Extraction)**: 최종 보고서에 파일명 대신 FASTA 헤더의 고유 식별자(e.g., accession number)를 샘플 ID로 사용하여 보고서의 가독성과 명확성을 향상시킵니다.
- **포괄적인 유전자 탐지**: 다음을 포함한 광범위한 중요 유전자를 식별합니다:
    - 항생제 내성 (AMR, Antimicrobial Resistance) 유전자
    - 플라스미드 복제원 (Plasmid Replicons)
    - 이동성 유전 인자 (MGEs, Mobile Genetic Elements)
- **고성능 비동기(Asynchronous) 처리**: Python의 `asyncio` 라이브러리를 활용하여 여러 BLAST 분석을 동시에 실행함으로써 총 분석 시간을 크게 단축합니다.
- **유연한 확장 구조 (Chain of Responsibility)**: "Chain of Responsibility" 디자인 패턴을 채택하여 분석 로직을 모듈화했습니다. 이 구조 덕분에 새로운 분석 유형을 파이프라인의 핵심 로직 수정 없이 깨끗하게 추가할 수 있습니다.
- **상세 로깅**: 전용 로거(logger)가 각 분석 모듈의 단계별 출력을 캡처하여 디버깅을 용이하게 하고 투명한 감사 추적(audit trail)을 제공합니다. `--verbose` 플래그는 실시간 콘솔 피드백도 제공합니다.
- **개발자 친화적인 코드베이스**: 소스 코드는 각 함수의 목적과 단계별 로직을 설명하는 구조화된 주석과 문서 문자열(docstring)로 광범위하게 문서화되어 있어, 새로운 기여자가 코드를 이해하고 수정하기 용이합니다.

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
│   │   ├── manager.py      #    - 전체 분석 파이프라인 조율 및 핸들러 체인 구성
│   │   ├── handler/        #    - "Chain of Responsibility" 패턴을 구현한 분석 핸들러 패키지
│   │   │   ├── __init__.py #    - 핸들러 모듈 초기화
│   │   │   ├── base.py     #    - 추상 핸들러 클래스 및 컨텍스트
│   │   │   ├── mlst.py     #    - MLST 분석 핸들러
│   │   │   ├── amr.py      #    - 항생제 내성(AMR) 분석 핸들러
│   │   │   ├── standard.py #    - 표준 BLAST 분석 핸들러
│   │   │   └── pathogen_finder.py # - PathogenFinder2 분석 핸들러
│   │   ├── blast_runner.py #    - BLAST 커맨드를 비동기적으로 실행
│   │   ├── pathogen_runner.py #    - PathogenFinder2 커맨드를 비동기적으로 실행
│   │   └── utils.py        #    - 헬퍼 함수 (e.g., 의존성 체크, 종 식별)
│   └── reporting/          #    - 리포트 생성 모듈
│       └── reporter.py     #    - 결과를 최종 텍스트 리포트로 컴파일
│
├── database/               # 2. 모든 분석 데이터베이스의 루트 디렉토리
│   ├── MLST_DB/            #    - 종별 MLST 데이터베이스
│   ├── resfinder_db/       #    - 항생제 내성 유전자 데이터베이스
│   ├── plasmidfinder_db/   #    - 플라스미드 복제원 데이터베이스
│   ├── mefinder_db/        #    - 이동성 유전 인자 데이터베이스
│   └── Pathogenfinder/     #    - PathogenFinder2 데이터베이스
│       ├── configs/        #    - PathogenFinder2 설정 파일
│       └── PathogenFinder2_dataset/ #    - PathogenFinder2 학습 데이터셋
│
├── genome/                 # 3. 입력 유전체 파일 위치
│   └── ...
├── analysis_results/       # 4. 리포트 및 결과물의 기본 출력 디렉토리
│
└── logs/                   # 5. 디버깅 및 감사 추적을 위한 상세 로그
```

---

## 4. 출력 구조 (Output Structure)

파이프라인은 `analysis_results`, `logs`, `blast_db_output` 디렉토리에 결과를 생성합니다. 모든 출력은 입력 유전체의 ID와 종(species)에 따라 체계적으로 구성됩니다.

예를 들어, `genome/test_id/klebsiella/A0018KP0093.fasta`를 입력으로 사용할 경우, 출력은 다음과 같이 생성됩니다:

```
analysis_results/
└── test_id/
    └── klebsiella/
        ├── Antimicrobial_Resistance/
        │   ├── blast_results.tsv
        │   └── combined_query.fasta
        ├── Mobile_Genetic_Elements/
        │   ├── blast_results.tsv
        │   └── combined_query.fasta
        ├── MLST/
        │   └── mlst_results.json  # <-- MLST 결과 (JSON 형식)
        ├── Plasmid_Replicons/
        │   ├── blast_results.tsv
        │   └── combined_query.fasta
        └── Final_ME_Report.txt      # <-- 최종 요약 보고서

logs/
└── test_id/
    └── klebsiella/
        ├── 2025-10-16_Pipeline_1_Pre-flight_Checks_1.log
        ├── 2025-10-16_MLST_1_Start_MLST_Workflow_1.log
        └── ... (각 단계별 상세 로그)

blast_db_output/
└── test_id/
    └── klebsiella/
        ├── A0018KP0093.nhr
        ├── A0018KP0093.nin
        └── A0018KP0093.nsq
```

### 주요 출력 파일 설명

-   **`analysis_results/{id}/{species}/Final_ME_Report.txt`**: 모든 분석 결과를 요약한 최종 텍스트 보고서입니다.
-   **`analysis_results/{id}/{species}/{analysis_name}/blast_results.tsv`**: 각 표준 분석에 대한 원시 BLAST 결과를 담고 있는 TSV 파일입니다.
-   **`analysis_results/{id}/{species}/MLST/mlst_results.json`**: MLST 분석 결과를 담고 있는 JSON 파일입니다. `tseemann/mlst`와 유사한 형식으로 ST, scheme, allele 프로파일을 포함합니다.
-   **`logs/{id}/{species}/`**: 파이프라인의 모든 단계에 대한 상세 로그 파일이 저장되는 디렉토리입니다. 디버깅에 유용합니다.
-   **`blast_db_output/{id}/{species}/`**: 입력 유전체로부터 생성된 BLAST 데이터베이스 파일들이 저장됩니다.

---

## 5. 코드 로직 및 모듈 설명

파이프라인의 핵심 로직은 "Chain of Responsibility" 디자인 패턴을 중심으로 구축되어, 각 모듈이 명확한 책임을 갖습니다.

- **`analysis/manager.py` (The Orchestrator)**
  - `AnalysisManager` 클래스는 전체 워크플로우를 총괄하는 오케스트레이터입니다.
  - `run_pipeline` 메소드는 파이프라인의 모든 단계를 순서대로 실행합니다: (1) 환경 설정, (2) BLAST DB 생성, (3) 분석 실행, (4) 리포팅, (5) 임시 파일 정리.
  - 가장 중요한 역할은 `analysis/handler/`에 정의된 **핸들러(handler)들을 조립하여 "책임 연쇄(Chain of Responsibility)"를 구성**하는 것입니다. `MLSTHandler`, `AMRHandler`, `PathogenFinder2Handler`, `StandardAnalysisHandler`를 순서대로 연결하여 체인을 만듭니다. 그 후, `config.py`에 명시된 모든 분석 요청을 체인의 첫 번째 핸들러에게 전달합니다.

- **`analysis/handler/` (The Handler Package)**
  - 이 패키지는 "Chain of Responsibility" 디자인 패턴을 구현하여 파이프라인의 유연성을 책임지는 핵심부입니다. 각 분석 유형은 별도의 핸들러 모듈로 분리되어 있습니다.
  - **동작 원리**: 분석 요청이 들어오면, 체인의 첫 번째 핸들러(e.g., `MLSTHandler`)가 요청을 확인합니다. 자신이 처리할 수 있는 요청이면 분석을 수행하고, 그렇지 않으면 다음 핸들러에게 요청을 그대로 전달합니다.
  - **주요 구성 요소**:
    - `handler/base.py`: 모든 핸들러가 공유하는 `AnalysisContext` 데이터 클래스와, 모든 핸들러가 상속해야 하는 추상 베이스 클래스 `AnalysisHandler`를 정의합니다.
    - `handler/mlst.py`: 복잡한 다단계 워크플로우를 가진 "특별(special)" MLST 분석을 처리하는 `MLSTHandler`를 포함합니다.
    - `handler/amr.py`: ABRicate와 유사하게 BLAST 결과를 파싱하고 요약하여 항생제 내성(AMR) 유전자를 식별하는 `AMRHandler`를 포함합니다.
    - `handler/pathogen_finder.py`: PathogenFinder2 분석을 위한 `PathogenFinder2Handler`를 포함합니다.
    - `handler/standard.py`: 체인의 마지막에 위치하며, 간단한 단일 BLAST 검색으로 처리할 수 있는 모든 "표준(standard)" 분석(e.g., Plasmid Replicons, MGEs)을 담당하는 `StandardAnalysisHandler`를 포함합니다.

- **`analysis/blast_runner.py` (The Asynchronous Worker)**
  - `makeblastdb`, `blastn`과 같은 NCBI BLAST+ 커맨드 라인 툴에 대한 비동기 래퍼(wrapper)를 제공합니다.
  - `asyncio.create_subprocess_exec`를 사용하여 여러 BLAST 작업을 병렬로 처리함으로써 파이프라인의 성능을 극대화합니다.

- **`config.py` (Central Configuration)**
  - `ANALYSES_TO_RUN` 딕셔너리는 파이프라인이 수행할 분석 목록을 정의하는 **컨트롤 센터**입니다. `manager`는 이 딕셔너리를 읽어 각 분석을 핸들러 체인에 전달합니다.

- **`utils.py`, `main.py`, `reporting.py`, `logger.py`**
  - 이 모듈들은 각각 의존성 확인 및 파라미터 준비(`utils`), CLI 인터페이스(`main`), 최종 리포트 생성(`reporting`), 상세 로그 기록(`logger`) 등 명확히 분리된 보조 기능을 수행합니다.

---

## 6. 시스템 요구사항 (System Requirements)

- **Python**: 3.9 이상
- **NCBI BLAST+ Suite**: `blastn`, `makeblastdb`, `blastdbcmd`가 설치되어 시스템의 PATH에서 접근 가능해야 합니다.
- **PathogenFinder2 Dependencies**:
  - `prodigal`: 유전체 예측 도구
  - `protT5`: 단백질 임베딩 도구
  - `diamond`: 고성능 단백질 비교 도구
- **Python Libraries**:
  - `pandas`
  - `biopython`

---

## 7. 설치 및 설정 (Installation and Setup)

**1. Repository 클론**
```bash
git clone <repository_url>
cd genome_analyzer
```

**2. Python 라이브러리 설치**
```bash
pip install pandas biopython
```

**3. 의존성 설치 확인**
파이프라인 실행 전에 모든 의존성이 설치되어 있는지 확인합니다:
```bash
# BLAST+ 도구 확인
which blastn makeblastdb blastdbcmd

# PathogenFinder2 도구 확인
which prodigal protT5 diamond
```

**4. 데이터베이스 설정**
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
├── mefinder_db/            # MGE 데이터베이스 파일 (*.fasta)
│   └── MGEdb_cds.fasta
│
└── Pathogenfinder/         # PathogenFinder2 데이터베이스
    ├── configs/            # PathogenFinder2 설정 파일
    │   ├── config_empty.json
    │   └── config_train.json
    └── PathogenFinder2_dataset/ # PathogenFinder2 학습 데이터셋
        └── ...
```

**5. PathogenFinder2 데이터베이스 설정**
PathogenFinder2를 사용하려면 다음과 같이 데이터베이스를 설정해야 합니다:

```bash
# PathogenFinder2 데이터베이스 다운로드
cd database/Pathogenfinder
# PathogenFinder2 공식 저장소에서 최신 데이터셋을 다운로드합니다
# (https://github.com/genomicepidemiology/PathogenFinder2)

# 데이터베이스 구조 확인
ls -la PathogenFinder2_dataset/
```

**6. 의존성 검증**
설치가 완료되었는지 확인합니다:
```bash
python -c "from src.analysis.utils import check_dependencies; check_dependencies()"
```

---

## 8. PathogenFinder2 통합 (PathogenFinder2 Integration)

PathogenFinder2는 병원체 유전체의 종 분류와 위험도 평가를 위한 고급 분석 도구입니다. 이 통합을 통해 파이프라인은 병원체 식별, 특성화, 위험도 평가를 자동화할 수 있습니다.

### 작동 원리

PathogenFinder2는 다음과 같은 단계로 작동합니다:

1. **의존성 확인**: `prodigal`, `protT5`, `diamond` 도구의 설치 여부를 확인합니다
2. **환경 설정**: 입력 유전체 경로와 출력 디렉토리를 포함한 설정 파일을 생성합니다
3. **실행**: PathogenFinder2를 실행하여 병원체 분석을 수행합니다
4. **결과 검증**: 생성된 결과 파일의 유효성을 검증하고 파싱합니다
5. **정리**: 임시 파일을 정리합니다

### 출력 형식 및 해석

PathogenFinder2는 다음과 같은 출력 파일을 생성합니다:

- **`pathogenfinder_results.tsv`**: 상세한 분석 결과 (TSV 형식)
- **`pathogenfinder_summary.txt`**: 요약 정보 및 위험도 평가

주요 결과 항목:
- **Pathogen Identification**: 병원체 종 식별 결과
- **Risk Assessment**: 위험도 등급 (High, Medium, Low)
- **Confidence Score**: 분신도 점수 (0-1)
- **Supporting Evidence**: 분석 근거 정보

### 문제 해결 (Troubleshooting)

**일반적인 문제 및 해결 방법:**

1. **의존성 오류**:
   ```
   Error: Missing PathogenFinder2 dependencies: prodigal, protT5, diamond
   ```
   **해결**: 모든 의존성이 설치되어 있는지 확인하고 PATH에 등록합니다

2. **데이터베이스 오류**:
   ```
   Error: PathogenFinder2 database not found
   ```
   **해결**: `database/Pathogenfinder/PathogenFinder2_dataset/` 디렉토리에 데이터베이스가 있는지 확인합니다

3. **메모리 부족**:
   ```
   Error: Out of memory during PathogenFinder2 execution
   ```
   **해결**: `-t` 옵션으로 스레드 수를 줄여 실행합니다 (예: `-t 2`)

4. **결과 파일 누락**:
   ```
   Error: Missing expected output files
   ```
   **해결**: 실행이 정상적으로 완료되었는지 확인하고, 출력 디렉토리 권한을 확인합니다

**디버깅 팁:**
- `--verbose` 플래그를 사용하여 상세한 로그를 확인합니다
- `logs/` 디렉토리에서 PathogenFinder2 관련 로그를 분석합니다
- 임시 디렉토리에 생성된 설정 파일을 검증합니다

---

## 9. 확장 가이드: 표준 분석 추가 (Expansion Guide: Standard Analysis)

단일 BLAST 검색으로 구성된 "표준" 분석(e.g., 새로운 Virulence Factor 데이터베이스)을 추가하는 것은 `config.py` 파일을 수정하는 것만으로 매우 간단하게 완료할 수 있습니다.

**1단계: 데이터베이스 추가**
새로운 FASTA 데이터베이스 파일을 `database/` 디렉토리 안의 전용 폴더에 배치합니다.

**2단계: `config.py` 설정 업데이트**
`src/config.py`를 열고 `ANALYSES_TO_RUN` 딕셔너리에 새 항목을 추가합니다.

```python
# in src/config.py
ANALYSES_TO_RUN = {
    "MLST_DB": "MLST",
    "virulence_db": "Virulence_Factors", # <-- 여기에 새 분석을 추가
    ...
}
```
이제 모든 준비가 끝났습니다. 파이프라인을 실행하면, 체인에 있는 다른 특별 핸들러가 요청을 처리하지 않는 한, `StandardAnalysisHandler`가 이 새로운 요청을 자동으로 감지하여 처리합니다.

---

## 10. 확장 가이드: 특별 분석 추가 (Expansion Guide: Special Analysis)

"Chain of Responsibility" 패턴 덕분에, MLST, AMR, PathogenFinder2와 같이 여러 단계의 복잡한 로직을 갖는 "특별" 분석을 추가하는 과정 또한 매우 체계적입니다.

**특별 핸들러 예시 (AMR, PathogenFinder2):**

`AMRHandler`나 `PathogenFinder2Handler`와 같은 특별 핸들러는 `StandardAnalysisHandler`가 처리하는 표준 분석과 다르게 다음과 같은 특징을 가집니다:

- **다단계 워크플로우**: 단순 BLAST 실행을 넘어, 결과를 파싱(`AMRHandler`)하거나 외부 도구 실행(`PathogenFinder2Handler`)과 같은 복합적인 단계를 포함합니다.
- **복잡한 의존성 관리**: `pandas`(`AMRHandler`)나 `prodigal`, `diamond`(`PathogenFinder2Handler`) 등 추가적인 라이브러리나 외부 도구를 관리합니다.
- **특화된 출력 생성**: 원시 BLAST 결과 외에, `amr_summary.json`(`AMRHandler`)과 같은 요약 파일을 생성합니다.
- **에러 처리**: 각 단계에서 발생할 수 있는 오류를 보다 체계적으로 처리합니다.

**1단계: `config.py`에 분석 정의**

표준 분석과 마찬가지로 `src/config.py`의 `ANALYSES_TO_RUN` 딕셔너리에 새로운 분석을 추가합니다.

```python
# in src/config.py
ANALYSES_TO_RUN = {
    "MLST_DB": "MLST",
    "Pathogenfinder": "Pathogen_Finder2",
    "resfinder_db": "Antimicrobial_Resistance",
    "serotype_db": "Serotyping", # <-- 새로운 특별 분석 추가
    ...
}
```

**2단계: `analysis/handler/` 패키지에 전용 핸들러 모듈 생성**

`src/analysis/handler/` 패키지 안에, 새로운 특별 분석을 위한 전용 핸들러 모듈(e.g., `serotyping.py`)을 작성합니다. 이 모듈의 클래스는 `AnalysisHandler`를 상속받아야 하며, `handle` 메소드 내에서 자신의 분석 이름(예: "Serotyping")을 확인하고, 일치할 경우 커스텀 워크플로우를 실행해야 합니다. `handler/mlst.py`의 `MLSTHandler`나 `handler/pathogen_finder.py`의 `PathogenFinder2Handler`를 템플릿으로 사용할 수 있습니다.

```python
# in src/analysis/handler/serotyping.py

from .base import AnalysisHandler
import asyncio

class SerotypingHandler(AnalysisHandler):
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        if analysis_name == "Serotyping":
            return asyncio.create_task(self._run_serotyping_workflow(params))
        else:
            return await super().handle(analysis_name, db_folder, params)

    async def _run_serotyping_workflow(self, params: dict):
        # Serotyping 워크플로우 구현
        ...
```

**3단계: `AnalysisManager`에서 핸들러 체인에 연결**

마지막으로, `src/analysis/manager.py`의 `run_pipeline` 메소드에서 새로 만든 핸들러를 체인에 연결합니다. `StandardAnalysisHandler`는 항상 체인의 가장 마지막에 위치해야 합니다.

```python
# in run_pipeline() method of AnalysisManager

            # Build the chain of responsibility
            standard_handler = StandardAnalysisHandler(context)
            pathogen_handler = PathogenFinder2Handler(context)
            amr_handler = AMRHandler(context)
            analysis_chain = MLSTHandler(context)

            # 체인 연결: MLST -> AMR -> PathogenFinder2 -> Standard
            analysis_chain.set_next(amr_handler).set_next(pathogen_handler).set_next(standard_handler)
```

**특별 분석 vs 표준 분석:**

| 특징 | 표준 분석 (StandardAnalysisHandler) | 특별 분석 (MLST, AMR, PathogenFinder2 등) |
|------|-----------------------------------|-----------------------------------|
| **워크플로우** | 단일 BLAST 검색 | 다단계 복합 워크플로우 (e.g., BLAST + 결과 파싱, 외부 도구 실행) |
| **의존성** | BLAST+만 필요 | BLAST+ 외 추가 라이브러리(e.g., pandas) 또는 외부 도구 필요 |
| **출력 처리** | 원시 BLAST 결과(TSV) 저장 | JSON 요약, 커스텀 리포트 등 추가적인 결과물 생성 |
| **에러 처리** | 기본 에러 처리 | 각 단계에 맞는 상세한 에러 처리 및 복구 로직 |
| **확장성** | `config.py` 수정만으로 가능 | 전용 핸들러 클래스 구현 필요 |

이 3단계를 통해, 파이프라인의 핵심 로직을 수정하지 않고도 새로운 특별 분석을 깨끗하고 모듈화된 방식으로 추가할 수 있습니다.

---

## 11. 설정 예시 (Configuration Examples)

**PathogenFinder2 설정 예시:**

**1. `ANALYSES_TO_RUN` 설정:**

```python
# in src/config.py
ANALYSES_TO_RUN = {
    # 특별 분석 (전용 핸들러 존재)
    "MLST_DB": "MLST",
    "Pathogenfinder": "Pathogen_Finder2",  # PathogenFinder2 활성화
    "resfinder_db": "Antimicrobial_Resistance",

    # 표준 분석 (StandardAnalysisHandler가 처리)
    "plasmidfinder_db": "Plasmid_Replicons",
    "mefinder_db": "Mobile_Genetic_Elements",
}
```

**2. 데이터베이스 설정 예시:**

```bash
# 데이터베이스 디렉토리 구조
database/
├── Pathogenfinder/
│   ├── configs/
│   │   ├── config_empty.json    # 기본 설정 파일
│   │   └── config_train.json    # 학습용 설정 파일
│   └── PathogenFinder2_dataset/ # 학습 데이터셋
│       ├── bacteria.faa         # 박테리아 단백질 데이터
│       ├── virus.faa            # 바이러스 단백질 데이터
│       └── fungi.faa            # 곰팡이 단백질 데이터
```

**3. 실행 예시:**

```bash
# PathogenFinder2를 포함한 전체 분석 실행
python src/main.py -g genome/sample.fasta -v

# 특정 분석만 실행 (PathogenFinder2 제외)
python src/main.py -g genome/sample.fasta -v --analyses MLST,Antimicrobial_Resistance

# PathogenFinder2만 실행
python src/main.py -g genome/sample.fasta -v --analyses Pathogen_Finder2
```

**4. 출력 디렉토리 구조:**

```bash
analysis_results/
├── MLST/
│   └── ... (MLST 결과)
├── Pathogen_Finder2/
│   ├── config.json              # PathogenFinder2 설정 파일
│   ├── pathogenfinder_results.tsv  # 상세 결과
│   ├── pathogenfinder_summary.txt  # 요약 정보
│   └── logs/                    # 분석 로그
├── Antimicrobial_Resistance/
│   └── ... (항생제 내성 결과)
└── ...
```

**5. 설정 파일 예시 (`config.json`):**

```json
{
  "input_genome": "/path/to/genome/sample.fasta",
  "output_dir": "/path/to/analysis_results/Pathogen_Finder2",
  "database_dir": "/path/to/database/Pathogenfinder",
  "threads": 4,
  "evalue_threshold": 1e-5
}
```

이 설정 예시들은 PathogenFinder2 통합을 위한 기본적인 구성을 보여줍니다. 실제 사용 시에는 데이터베이스 구조와 환경에 맞게 설정을 조정해야 합니다.
