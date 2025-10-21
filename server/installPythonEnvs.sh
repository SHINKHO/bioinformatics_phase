#!/bin/bash

# 이 스크립트는 Miniconda(Conda)와 pyenv를 사용하여 Python 개발 환경의 '기반'을 구축합니다.
# Conda의 기본 채널을 'conda-forge'로 설정하고 'defaults' 채널은 제거합니다.
set -e

echo "Python & Conda 기본 환경 설치를 시작합니다."
echo "========================================================"

# --- Conda 설치 및 설정 섹션 ---
if command -v conda &> /dev/null; then
    echo "Conda가 이미 설치되어 있습니다. 설치를 건너뜁니다."
else
    echo "[1/2] Miniconda (Conda)를 설치합니다..."
    curl -fsSL -o ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash ~/miniconda.sh -b -u -p "$HOME/miniconda3"
    rm ~/miniconda.sh
    "$HOME/miniconda3/bin/conda" init bash
    
    echo "[2/2] Conda 채널을 설정합니다 ('conda-forge'를 기본으로, 'defaults' 제거)..."
    "$HOME/miniconda3/bin/conda" config --add channels conda-forge
    "$HOME/miniconda3/bin/conda" config --set channel_priority strict
    "$HOME/miniconda3/bin/conda" config --remove channels defaults
fi

echo "========================================================"
echo "✅ Conda 기본 환경 설정이 완료되었습니다."
echo "❗️ 중요: 변경 사항을 적용하려면 셸을 다시 시작해야 합니다."
echo "   오케스트레이터가 이 과정을 자동으로 처리합니다."
echo "========================================================"
