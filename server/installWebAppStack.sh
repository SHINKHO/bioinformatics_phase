#!/bin/bash

# 이 스크립트는 Conda를 사용하여 Django 및 uWSGI를 포함하는 특정 환경을 설정합니다.
# 이 스크립트는 오케스트레이터에 의해 호출되며, 비대화형으로 실행됩니다.
set -e

# --- 설정 변수 (init.sh로부터 환경 변수로 전달받음) ---
if [ -z "$ENV_NAME" ]; then
    echo "오류: ENV_NAME 환경 변수가 설정되지 않았습니다."
    echo "이 스크립트는 init.sh에 의해 호출되어야 합니다."
    exit 1
fi

PYTHON_VERSION="3.11"
PACKAGES_TO_INSTALL="django uwsgi"
# -----------------------------

echo "Django/uWSGI 웹 앱 스택 설정을 시작합니다."
echo "대상 환경: $ENV_NAME"
echo "========================================================"

# Conda 명령어 경로를 명시적으로 지정
CONDA_CMD="$HOME/miniconda3/bin/conda"

if ! command -v "$CONDA_CMD" &> /dev/null; then
    echo "오류: Conda 명령어를 찾을 수 없습니다. Python 기본 환경이 올바르게 설치되지 않았을 수 있습니다."
    exit 1
fi

# 1단계: Conda 환경 생성 (존재하지 않는 경우)
if ! "$CONDA_CMD" env list | grep -q "^${ENV_NAME}\s"; then
    echo "[1/2] Python v${PYTHON_VERSION}으로 '$ENV_NAME' Conda 환경을 생성합니다..."
    "$CONDA_CMD" create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
else
    echo "[1/2] '$ENV_NAME' 환경이 이미 존재합니다. 환경 생성을 건너뜁니다."
fi

# 2단계: Django 및 uWSGI 설치 (항상 실행)
# 환경이 이미 존재하더라도, 패키지가 누락되었을 수 있으므로 이 단계는 항상 실행합니다.
# pip install은 멱등성이 있어 이미 설치된 패키지는 자동으로 건너뜁니다.
echo "[2/2] '$ENV_NAME' 환경에 Django와 uWSGI를 설치/업데이트합니다..."
"$CONDA_CMD" run -n "$ENV_NAME" pip install $PACKAGES_TO_INSTALL

echo "========================================================"
echo "✅ '$ENV_NAME' 환경에 웹 앱 스택 설치/업데이트가 완료되었습니다!"
echo "   다음 명령어로 환경을 활성화할 수 있습니다: conda activate $ENV_NAME"
echo "========================================================"