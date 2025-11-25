# 단계별 디버그 정보를 파일로 저장하는 로거.
import logging
from pathlib import Path
from datetime import datetime

class Logger:
    # 단계별 디버그 정보를 개별 파일에 저장하는 간단한 로거.

    def __init__(self, log_dir: Path):
        # 로거 초기화. In: log_dir(Path) / Out: None
        self.log_dir = log_dir
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.log_counts = {}

    def log_step(self, analysis_type: str, step_name: str, content: str, extension: str = "log"):
        # 주어진 내용을 구조화된 이름의 파일에 로그로 기록. In: analysis_type(str), step_name(str), content(str), extension(str) / Out: None
        try:
            safe_step_name = "".join(c for c in step_name if c.isalnum() or c in ('_', '-')).rstrip()
            
            log_key = (analysis_type, safe_step_name)
            if log_key not in self.log_counts:
                self.log_counts[log_key] = 0
            self.log_counts[log_key] += 1
            count = self.log_counts[log_key]
            
            date_str = datetime.now().strftime("%Y-%m-%d")
            log_file_name = f"{date_str}_{analysis_type}_{safe_step_name}_{count}.{extension}"
            log_file = self.log_dir / log_file_name

            with open(log_file, "w", encoding="utf-8") as f:
                f.write(content)
        except Exception as e:
            print(f"Failed to write log for step '{step_name}'. Error: {e}")


def setup_run_logger(run_log_dir: Path):
    # (미사용) 콘솔 출력을 위한 메인 실행 로거 설정. In: run_log_dir(Path) / Out: logging.Logger
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.INFO)
    
    if not logger.handlers:
        c_handler = logging.StreamHandler()
        c_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        c_handler.setFormatter(c_format)
        logger.addHandler(c_handler)

    return logger
