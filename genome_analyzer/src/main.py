# 게놈 분석 파이프라인 메인 실행 스크립트.
import argparse
import asyncio
from pathlib import Path

# Project-level module imports
from config import DEFAULT_RESULTS_DIR
from analysis.manager import AnalysisManager

def main():
    # 커맨드라인 인자를 파싱하고 분석 파이프라인을 시작. In: None / Out: None
    parser = argparse.ArgumentParser(
        description="A comprehensive genome analysis pipeline for molecular epidemiology.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "genome_file",
        type=Path,
        help="Path to the input genome file in FASTA format."
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help=f"Path to the output directory.\n(default: {DEFAULT_RESULTS_DIR})"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose mode to see detailed progress logs."
    )
    
    args = parser.parse_args()
    
    manager = AnalysisManager(
        genome_file=args.genome_file,
        results_dir=args.output,
        verbose=args.verbose
    )
    
    asyncio.run(manager.run_pipeline())

if __name__ == "__main__":
    main()