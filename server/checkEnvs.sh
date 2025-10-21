#!/bin/bash

# 이 스크립트는 오케스트레이터에 의해 실행되며, 시스템 환경을 점검하고
# 결과를 변수로 설정합니다. 화면 출력은 최소화합니다.

# --- 상태 변수 초기화 ---
NGINX_MISSING=false
PYTHON_ENVS_MISSING=false
WEBAPP_STACK_MISSING=false

# --- 점검 1: NGINX ---
if ! command -v nginx &> /dev/null; then
    NGINX_MISSING=true
fi

# --- 점검 2: Python/Conda 기본 환경 ---
# Conda 명령어 경로를 다른 스크립트와 일치시킵니다.
CONDA_CMD="$HOME/miniconda3/bin/conda"

# conda 명령어가 없으면 Python 환경 전체가 없다고 간주합니다.
if ! command -v "$CONDA_CMD" &> /dev/null; then
    PYTHON_ENVS_MISSING=true
    WEBAPP_STACK_MISSING=true # Conda가 없으면 당연히 웹 앱 스택도 없습니다.
    return # 더 이상 점검할 필요가 없으므로 종료합니다.
fi

# --- 점검 3: Django/uWSGI Conda 환경 ---
# $ENV_NAME 변수는 init.sh에서 source될 때 전달받아야 합니다.
if [ -z "$ENV_NAME" ]; then
    echo "checkEnvs.sh: \$ENV_NAME 변수가 설정되지 않았습니다. 점검을 건너뜁니다." >&2
    WEBAPP_STACK_MISSING=true
    return 
fi

# 3-1: 환경 자체가 존재하는지 확인
if ! "$CONDA_CMD" env list | grep -q "^${ENV_NAME}\s"; then
    WEBAPP_STACK_MISSING=true
else
    # 3-2: 환경이 존재하면, 패키지(Django, uWSGI)가 설치되었는지 확인
    # 'pip list' 결과를 변수에 저장하여 한 번만 실행합니다.
    # 2>/dev/null은 conda run 실행 시 발생하는 stderr 출력을 숨깁니다.
    INSTALLED_LIST=$("$CONDA_CMD" run -n "$ENV_NAME" pip list 2>/dev/null)

    # 'pip install django'는 'Django'를, 'pip install uwsgi'는 'uWSGI'를 설치합니다.
    # 하나라도 없으면 WEBAPP_STACK_MISSING=true로 설정합니다.
    if ! echo "$INSTALLED_LIST" | grep -q -E "^Django\s"; then
        WEBAPP_STACK_MISSING=true
    elif ! echo "$INSTALLED_LIST" | grep -q -E "^uWSGI\s"; then
        WEBAPP_STACK_MISSING=true
    fi
fi