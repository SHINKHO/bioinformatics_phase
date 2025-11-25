
import asyncio
from abc import ABC, abstractmethod
from pathlib import Path
from dataclasses import dataclass

from logger import Logger

# --- 데이터 컨텍스트 ---

@dataclass
class AnalysisContext:
    # 모든 핸들러가 공유하는 데이터, 도구, 설정을 담는 데이터 클래스.
    genome_db_path: Path
    results_dir: Path
    temp_dir: Path
    logger: Logger
    verbose: bool
    results_data: dict
    genome_id: str
    species: str

# --- 핸들러 ABC ---

class AnalysisHandler(ABC):
    # 모든 분석 핸들러의 추상 기본 클래스 (ABC).
    def __init__(self, context: AnalysisContext):
        self._next_handler: AnalysisHandler | None = None
        self._context = context

    def set_next(self, handler: 'AnalysisHandler') -> 'AnalysisHandler':
        # 핸들러를 체인의 다음 핸들러와 연결. In: handler(AnalysisHandler) / Out: handler(AnalysisHandler)
        self._next_handler = handler
        return handler

    @abstractmethod
    async def handle(self, analysis_name: str, db_folder: str, params: dict) -> asyncio.Task | None:
        # 분석 요청 처리. In: analysis_name(str), db_folder(str), params(dict) / Out: asyncio.Task 또는 None
        if self._next_handler:
            return await self._next_handler.handle(analysis_name, db_folder, params)
        return None
