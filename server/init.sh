#!/bin/bash

# 이 스크립트는 웹 애플리케이션 스택 설치 과정을 총괄합니다.
# 1. checkEnvs.sh를 실행하여 현재 시스템 상태를 점검합니다.
# 2. 점검 결과에 따라 필요한 설치 스크립트를 순서대로 호출합니다.

# 실행 중 오류가 발생하면 즉시 중단
set -e

# --- 스크립트 실행 경로 설정 ---
# 이 스크립트가 위치한 디렉토리를 기준으로 다른 스크립트를 찾습니다.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# --- 중앙 설정 변수 ---
# Conda 환경 이름을 여기서 중앙 관리합니다.
# export를 통해 이 스크립트가 실행하는 모든 하위 셸(bash)에 이 변수가 전달됩니다.
export ENV_NAME="blastn_test"

# --- 1단계: 초기 환경 점검 ---
echo "========================================================"
echo "1. 시스템 환경 점검을 시작합니다 (대상 환경: $ENV_NAME)..."
echo "========================================================"
# source는 현재 셸에서 실행되므로 $ENV_NAME 변수가 checkEnvs.sh에 공유됩니다.
source "${SCRIPT_DIR}/checkEnvs.sh"

# --- 2단계: 점검 결과에 따른 설치 진행 ---
echo ""
echo "========================================================"
echo "2. 점검 결과를 바탕으로 설치를 진행합니다..."
echo "========================================================"

# NGINX가 없는 경우
if [ "$NGINX_MISSING" = true ]; then
    echo "-> NGINX가 설치되지 않았습니다. 설치를 시작합니다..."
    bash "${SCRIPT_DIR}/installNginx.sh"
    echo "-> NGINX 설치 완료."
    echo ""
fi

# Python/Conda 기본 환경이 없는 경우
if [ "$PYTHON_ENVS_MISSING" = true ]; then
    echo "-> Python/Conda 기본 환경이 설정되지 않았습니다. 설치를 시작합니다..."
    bash "${SCRIPT_DIR}/installPythonEnvs.sh"
    
    echo ""
    echo "-> Python 환경 설치가 완료되었습니다. 변경 사항을 적용하기 위해 스크립트를 재시작합니다..."
    # 셸을 재시작하여 PATH 변경 등을 적용한 후, 이 스크립트를 다시 실행합니다.
    # 이렇게 해야 이후 단계에서 'conda' 명령어를 인식할 수 있습니다.
    exec "$SHELL" -c "cd '$SCRIPT_DIR' && bash '$0'"
fi

# 웹 앱 스택(Django/uWSGI) 설치
# 이 코드는 Python/Conda 설치가 완료된 후에만 도달합니다.
if [ "$WEBAPP_STACK_MISSING" = true ]; then
    echo "-> Django/uWSGI 환경이 없거나 패키지가 누락되었습니다. 설치/업데이트를 시작합니다..."
    bash "${SCRIPT_DIR}/installWebAppStack.sh"
    echo "-> Django/uWSGI 환경 설치 완료."
    echo ""
fi

# --- 최종 점검 ---
# 모든 설치가 끝난 후, 최종적으로 환경을 다시 점검합니다.
echo ""
echo "========================================================"
echo "3. 최종 점검을 수행합니다..."
echo "========================================================"
# 설치 시도 후 상태가 변경되었을 수 있으므로, checkEnvs를 다시 source하여 변수를 최신화합니다.
source "${SCRIPT_DIR}/checkEnvs.sh"

# --- 최종 결과 ---
echo ""
echo "========================================================"
echo "최종 설치 결과"
echo "========================================================"
if [ "$NGINX_MISSING" = false ] && [ "$PYTHON_ENVS_MISSING" = false ] && [ "$WEBAPP_STACK_MISSING" = false ]; then
    echo "✅ SUCCESS: 모든 필수 구성 요소가 성공적으로 설치 및 설정되었습니다."
else
    echo "❌ ERROR: 일부 구성 요소가 여전히 누락되었습니다. 로그를 확인해주세요."
    [ "$NGINX_MISSING" = true ] && echo "   - NGINX 누락"
    [ "$PYTHON_ENVS_MISSING" = true ] && echo "   - Python/Conda 기반 환경 누락"
    [ "$WEBAPP_STACK_MISSING" = true ] && echo "   - Django/uWSGI 환경 누락"
    exit 1
fi
echo "========================================================"