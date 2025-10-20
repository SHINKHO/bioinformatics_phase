# Genome Analyzer Pipeline

## 1. 소개 (Introduction)

Genome Analyzer Pipeline은 박테리아 유전체 서열의 분자 역학(Molecular Epidemiology) 분석을 자동화하기 위해 설계된 Python 기반의 데이터 중심(data-driven) 파이프라인입니다. `pipelines/` 폴더의 `*.json` 설정 파일을 통해 전체 워크플로우를 정의하며, FASTA 파일을 입력받아 MLST, 항생제 내성 유전자 탐지 등 일련의 생물정보학 분석을 수행하고, 그 결과를 하나의 포괄적인 요약 보고서로 통합합니다.

이 파이프라인의 주요 목표는 반복적인 분석 작업을 자동화하여 연구 효율성을 높이고, 박테리아 유전체에 대한 빠르고 재현 가능한(reproducible) 인사이트를 제공하는 것입니다.

---

## 2. 주요 기능 (Key Features)

- **데이터 중심 아키텍처 (Data-Driven Architecture)**: `pipelines/*.json` 파일을 통해 실행할 분석 단계, 각 단계의 입출력 및 파라미터를 명시적으로 정의합니다. 이를 통해 코드 수정 없이 파이프라인의 동작을 유연하게 변경할 수 있습니다.
- **자동 MLST 분석**: 7-gene MLST 스킴을 기반으로 housekeeping gene을 자동으로 식별하고, allele 번호를 결정하며, 최종 Sequence Type (ST)을 부여합니다.
- **PathogenFinder2 통합**: 병원체 식별 및 특성화를 위한 PathogenFinder2 분석을 지원하여 병원체 유전체의 종 분류와 위험도 평가를 수행합니다.
- **동적 종(Species) 식별**: 입력된 유전체의 FASTA 헤더에서 박테리아 종을 지능적으로 감지하여, 해당 종에 적합한 MLST 데이터베이스를 동적으로 선택하여 분석합니다.
- **지능적인 샘플 ID 추출 (Intelligent Sample ID Extraction)**: 최종 보고서에 파일명 대신 FASTA 헤더의 고유 식별자(e.g., accession number)를 샘플 ID로 사용하여 보고서의 가독성과 명확성을 향상시킵니다.
- **포괄적인 유전자 탐지**: 다음을 포함한 광범위한 중요 유전자를 식별합니다:
    - 항생제 내성 (AMR, Antimicrobial Resistance) 유전자
    - 플라스미드 복제원 (Plasmid Replicons)
    - 이동성 유전 인자 (MGEs, Mobile Genetic Elements)
- **고성능 비동기(Asynchronous) 처리**: Python의 `asyncio` 라이브러리를 활용하여 여러 분석 단계를 동시에 실행함으로써 총 분석 시간을 크게 단축합니다.
- **유연한 확장 구조 (Abstract Handler)**: 모든 분석 모듈은 `AnalysisHandler` 추상 클래스를 상속받아 구현됩니다. 이 구조 덕분에 새로운 분석 유형을 파이프라인의 핵심 로직 수정 없이 깨끗하고 독립적인 모듈로 추가할 수 있습니다.
- **상세 로깅**: 전용 로거(logger)가 각 분석 모듈의 단계별 출력을 캡처하여 디버깅을 용이하게 하고 투명한 감사 추적(audit trail)을 제공합니다. `--verbose` 플래그는 실시간 콘솔 피드백도 제공합니다.
- **개발자 친화적인 코드베이스**: 소스 코드는 각 함수의 목적과 단계별 로직을 설명하는 구조화된 주석과 문서 문자열(docstring)로 광범위하게 문서화되어 있어, 새로운 기여자가 코드를 이해하고 수정하기 용이합니다.

---

## 3. 프로젝트 구조 (Project Structure)

프로젝트는 소스 코드, 데이터베이스, 입출력 데이터를 위한 별도의 디렉토리로 구성되어 명확성과 유지보수성을 높입니다.

```
genome_analyzer/
├── pipelines/              # 0. 파이프라인 워크플로우 정의 폴더
│   └── default.json        #    - 기본 파이프라인
├── src/                    # 1. 메인 소스 코드
│   ├── main.py             #    - CLI 진입점, argument 파싱, pipeline 선택
│   ├── pipeline_service.py #    - pipeline.json을 읽어 전체 워크플로우를 조율하는 오케스트레이터
│   ├── config.py           #    - 경로 등 레거시 설정 (점차 pipeline 정의 파일로 통합)
│   ├── logger.py           #    - 단계별 디버그 로거
│   ├── analysis/           #    - 핵심 분석 로직 모듈
│   │   ├── context.py      #    - 파이프라인 단계 간 상태를 공유하는 컨텍스트
│   │   ├── handler/        #    - 분석 유형별 로직을 구현한 핸들러 패키지
│   │   │   ├── __init__.py #    - 핸들러 모듈 초기화
│   │   │   ├── base.py     #    - 모든 핸들러가 상속하는 추상 베이스 클래스
│   │   │   ├── mlst.py     #    - MLST 분석 핸들러
│   │   │   ├── amr.py      #    - 항생제 내성(AMR) 분석 핸들러
│   │   │   ├── standard.py #    - 표준 BLAST 분석 핸들러
│   │   │   └── pathogen_finder.py # - PathogenFinder2 분석 핸들러
│   │   ├── blast_runner.py #    - BLAST 커맨드를 비동기적으로 실행
│   │   ├── pathogen_runner.py # - PathogenFinder2 커맨드를 비동기적으로 실행
│   │   └── utils.py        #    - 헬퍼 함수 (e.g., 의존성 체크, 종 식별)
│   └── reporting/          #    - 리포트 생성 모듈
│       └── reporter.py     #    - 결과를 최종 텍스트 리포트로 컴파일
│
├── database/               # 2. 모든 분석 데이터베이스의 루트 디렉토리
│   └── ...
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
        │   └── mlst_results.json
        ├── Plasmid_Replicons/
        │   ├── blast_results.tsv
        │   └── combined_query.fasta
        └── Final_ME_Report.txt

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

-   `analysis_results/{id}/{species}/Final_ME_Report.txt`: 모든 분석 결과를 요약한 최종 텍스트 보고서입니다.
-   `analysis_results/{id}/{species}/{analysis_name}/blast_results.tsv`: 각 표준 분석에 대한 원시 BLAST 결과를 담고 있는 TSV 파일입니다.
-   `analysis_results/{id}/{species}/MLST/mlst_results.json`: MLST 분석 결과를 담고 있는 JSON 파일입니다. `tseemann/mlst`와 유사한 형식으로 ST, scheme, allele 프로파일을 포함합니다.
-   `logs/{id}/{species}/`: 파이프라인의 모든 단계에 대한 상세 로그 파일이 저장되는 디렉토리입니다. 디버깅에 유용합니다.
-   `blast_db_output/{id}/{species}/`: 입력 유전체로부터 생성된 BLAST 데이터베이스 파일들이 저장됩니다.

---

## 5. 코드 로직 및 모듈 설명

파이프라인의 핵심 로직은 `pipelines` 폴더 내의 JSON 파일에 정의된 워크플로우를 동적으로 실행하는 데이터 중심(data-driven) 아키텍처를 기반으로 합니다.

- **`pipelines/*.json` (The Workflow Definition)**
  - 파이프라인이 수행할 **분석 단계(step)의 목록**을 정의하는 JSON 파일입니다.
  - 각 단계는 사용할 핸들러 모듈, 입력/출력 경로, 그리고 필요한 파라미터를 명시합니다.
  - 이 파일을 수정하는 것만으로 코드 변경 없이 파이프라인의 순서와 구성을 변경할 수 있습니다.

- **`pipeline_service.py` (The Data-Driven Orchestrator)**
  - `PipelineService` 클래스는 전체 워크플로우를 총괄하는 오케스트레이터입니다.
  - `run_analysis` 메소드는 지정된 파이프라인 JSON 파일을 읽고, 정의된 각 단계를 순서대로 실행합니다.
  - 각 단계마다 `module`에 명시된 핸들러 클래스를 동적으로 임포트하고, 해당 단계의 설정(`step_config`)과 공유 컨텍스트(`WorkflowContext`)를 전달하여 핸들러 인스턴스를 생성한 후, `execute()` 메소드를 호출합니다.

- **`analysis/handler/` (The Handler Package)**
  - 이 패키지는 각 분석 유형의 실제 로직을 구현하는 독립적인 '플러그인' 모음입니다.
  - **동작 원리**: `PipelineService`가 `pipelines` 폴더의 JSON 파일에 명시된 핸들러를 직접 지정하여 실행하므로, 핸들러 간의 직접적인 연결(체인)은 없습니다. 각 핸들러는 자신에게 주어진 임무만 수행하면 됩니다.
  - **주요 구성 요소**:
    - `handler/base.py`: 모든 핸들러가 상속해야 하는 추상 베이스 클래스 `AnalysisHandler`를 정의합니다. 이 클래스는 `__init__`에서 입출력 경로를 자동으로 처리하고, 파일 읽기/쓰기를 위한 `_read_input()` / `_write_output()`과 같은 공통 유틸리티 메소드를 제공하여 핸들러 개발을 표준화하고 단순화합니다.
    - `handler/mlst.py`, `handler/amr.py`, `handler/pathogen_finder.py`: 각각 MLST, AMR, PathogenFinder2와 같이 복잡한 로직을 가진 "특별(special)" 분석을 처리합니다.
    - `handler/standard.py`: 간단한 단일 BLAST 검색으로 처리할 수 있는 모든 "표준(standard)" 분석(e.g., Plasmid Replicons, MGEs)을 담당합니다. `pipelines` 폴더의 JSON 파일에서 이 핸들러를 재사용하여 다양한 BLAST 기반 분석을 수행할 수 있습니다.

- **`analysis/blast_runner.py` (The Asynchronous Worker)**
  - `makeblastdb`, `blastn`과 같은 NCBI BLAST+ 커맨드 라인 툴에 대한 비동기 래퍼(wrapper)를 제공합니다.
  - `asyncio.create_subprocess_exec`를 사용하여 여러 BLAST 작업을 병렬로 처리함으로써 파이프라인의 성능을 극대화합니다.

- **`utils.py`, `main.py`, `reporting.py`, `logger.py`**
  - 이 모듈들은 각각 의존성 확인 및 파라미터 준비(`utils`), CLI 인터페이스(`main`), 최종 리포트 생성(`reporting`), 상세 로그 기록(`logger`) 등 명확히 분리된 보조 기능을 수행합니다.

---

## 6. 시스템 요구사항 (System Requirements)

- **Python**: 3.9+
- **NCBI BLAST+ Suite**: `blastn`, `makeblastdb`, `blastdbcmd`
- **PathogenFinder2 Dependencies**: `prodigal`, `protT5`, `diamond`
- **Python Libraries**: `pandas`, `biopython`

---

## 7. 설치 및 설정 (Installation and Setup)

**1. Repository 클론**
```bash
git clone https://github.com/your-username/genome_analyzer.git
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
    └── PathogenFinder2_dataset/ # PathogenFinder2 학습 데이터셋
```

**5. PathogenFinder2 데이터베이스 설정**
PathogenFinder2를 사용하려면 [PathogenFinder2 공식 저장소](https://github.com/genomicepidemiology/PathogenFinder2)에서 최신 데이터셋을 다운로드하여 `database/Pathogenfinder/PathogenFinder2_dataset`에 배치해야 합니다.

**6. 의존성 검증**
설치가 완료되었는지 확인합니다:
```bash
python -c "from src.analysis.utils import check_dependencies; check_dependencies()"
```

---

## 8. PathogenFinder2 통합 (PathogenFinder2 Integration)

파이프라인은 `PathogenFinder2Handler`를 통해 PathogenFinder2 분석을 통합합니다. 이 핸들러는 `pipelines` 설정 파일에서 `analysis.handler.pathogen_finder.PathogenFinder2Handler` 모듈로 지정될 때 트리거됩니다. 핸들러는 `pathogenfinder2` 커맨드 라인 도구를 적절한 파라미터와 함께 실행하고 결과를 캡처합니다. 또한, `prodigal`, `protT5`, `diamond`와 같은 필수 의존성이 시스템의 PATH에 있는지 확인하는 종속성 검사를 수행합니다.

---

## 9. 확장 가이드: 새로운 분석 추가하기

리팩토링된 아키텍처에서는 "표준" 분석과 "특별" 분석의 추가 방법이 거의 동일하며 매우 간단합니다. 모든 것은 `pipelines` 폴더 내의 JSON 파일을 수정하고, 필요에 따라 전용 핸들러를 만드는 두 단계로 귀결됩니다.

**1단계: (선택 사항) 전용 핸들러 모듈 생성**

- **단순 BLAST 분석**: 새로운 Virulence Factor 데이터베이스를 검색하는 것처럼 간단한 BLAST 분석을 추가하는 경우, 이 단계를 건너뛸 수 있습니다. 기존의 `analysis.handler.standard.StandardAnalysisHandler`를 재사용하면 됩니다.
- **복잡한 분석**: AMR이나 MLST처럼 다단계 워크플로우나 결과 파싱이 필요한 경우, `src/analysis/handler/` 디렉토리에 `my_new_analysis.py`와 같은 새 파일을 만듭니다. `AnalysisHandler`를 상속하는 클래스를 작성하고 `execute(self)` 메소드에 로직을 구현합니다. `amr.py`나 `mlst.py`를 훌륭한 템플릿으로 사용할 수 있습니다.

```python
# in src/analysis/handler/my_new_analysis.py
from .base import AnalysisHandler

class MyNewAnalysisHandler(AnalysisHandler):
    async def execute(self):
        self.logger.log_step(self.analysis_name, "Start", "Starting my new analysis.")
        # ... Your custom logic here ...
        # Use self.inputs, self.outputs, self.context
        # Use self._read_input() and self._write_output()
        self.logger.log_step(self.analysis_name, "End", "My new analysis complete.")
```

**2단계: `pipelines/*.json`에 새로운 단계(Step) 추가**

사용할 `pipelines` 폴더의 JSON 파일을 열고 `steps` 배열에 새로운 분석 단계를 위한 JSON 객체를 추가합니다.

```json
// in pipelines/my_pipeline.json
{
  "pipeline_name": "Standard Bacterial ME Analysis",
  "steps": [
    {
      "name": "MLST_Analysis",
      "enabled": true,
      "module": "analysis.handler.mlst.MLSTHandler",
      ...
    },
    {
      "name": "My_New_Analysis",
      "enabled": true,
      "module": "analysis.handler.my_new_analysis.MyNewAnalysisHandler", // 1단계에서 만든 핸들러
      "inputs": {
        "genome_db": "{context.genome_db_path}",
        "my_db": "database/my_new_db_folder" // 필요한 입력 파일
      },
      "outputs": {
        "results_file": "{context.results_dir}/My_New_Analysis/results.json" // 결과 파일
      },
      "parameters": {
        "my_param": "some_value" // 분석에 필요한 추가 파라미터
      }
    },
    // ... 다른 단계들
  ]
}
```

- **`module`**: 사용할 핸들러의 전체 Python 경로를 지정합니다. 단순 BLAST 분석의 경우 `analysis.handler.standard.StandardAnalysisHandler`를 지정하면 됩니다.
- **`inputs` / `outputs`**: 이 단계에서 사용할 입력 및 출력 파일의 경로를 지정합니다. `{context.*}` 플레이스홀더를 사용하여 동적으로 경로를 생성할 수 있습니다.
- **`parameters`**: 분석에 필요한 추가적인 설정값(e.g., BLAST의 identity threshold)을 전달합니다.

이 두 단계만으로 파이프라인의 핵심 로직(`pipeline_service.py`)을 전혀 수정하지 않고도 새로운 분석 기능을 깨끗하고 모듈화된 방식으로 추가할 수 있습니다.

---

## 10. 설정 예시 (Configuration Examples)

다음은 `pipelines/default.json` 파일의 예시입니다. 이 파일은 파이프라인의 각 단계를 정의합니다.

```json
{
  "pipeline_name": "Standard Bacterial ME Analysis",
  "steps": [
    {
      "name": "MLST_Analysis",
      "enabled": true,
      "module": "analysis.handler.mlst.MLSTHandler",
      "inputs": {
        "genome_file": "{context.initial_genome_file}",
        "species_db": "{context.database_root}/MLST_DB/{context.species}"
      },
      "outputs": {
        "results_json": "{context.results_dir}/MLST/mlst_results.json",
        "novel_alleles": "{context.results_dir}/MLST/novel_alleles.fasta"
      },
      "parameters": {}
    },
    {
      "name": "AMR_Analysis",
      "enabled": true,
      "module": "analysis.handler.amr.AMRHandler",
      "inputs": {
        "genome_db": "{context.genome_db_path}",
        "amr_db_folder": "resfinder_db"
      },
      "outputs": {
        "summary_json": "{context.results_dir}/Antimicrobial_Resistance/amr_summary.json"
      },
      "parameters": {
        "analysis_name": "Antimicrobial_Resistance"
      }
    },
    {
      "name": "Plasmid_Analysis",
      "enabled": true,
      "module": "analysis.handler.standard.StandardAnalysisHandler",
      "inputs": {
        "genome_db": "{context.genome_db_path}",
        "db_folder": "plasmidfinder_db"
      },
      "outputs": {
        "blast_results": "{context.results_dir}/Plasmid_Replicons/blast_results.tsv"
      },
      "parameters": {
        "analysis_name": "Plasmid_Replicons"
      }
    },
    {
      "name": "MGE_Analysis",
      "enabled": true,
      "module": "analysis.handler.standard.StandardAnalysisHandler",
      "inputs": {
        "genome_db": "{context.genome_db_path}",
        "db_folder": "mefinder_db"
      },
      "outputs": {
        "blast_results": "{context.results_dir}/Mobile_Genetic_Elements/blast_results.tsv"
      },
      "parameters": {
        "analysis_name": "Mobile_Genetic_Elements"
      }
    }
  ]
}
```
