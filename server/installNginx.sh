#!/bin/bash

# 이 스크립트는 NGINX 공식 문서(nginx.org)를 기반으로 Ubuntu 시스템에 NGINX를 설치합니다.
# 실행 중 오류가 발생하면 즉시 중단하고, 롤백 함수를 실행하도록 설정합니다.
set -e

# 롤백 함수: 스크립트 실패 시 생성된 파일들을 삭제하여 시스템을 원상 복구합니다.
rollback() {
  echo ""
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "오류가 발생하여 설치를 중단하고 롤백을 시작합니다..."
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  # NGINX 저장소 pinning 설정 파일 삭제
  if [ -f /etc/apt/preferences.d/99nginx ]; then
    echo "삭제 중: /etc/apt/preferences.d/99nginx"
    sudo rm -f /etc/apt/preferences.d/99nginx
  fi

  # NGINX 저장소 설정 파일 삭제
  if [ -f /etc/apt/sources.list.d/nginx.list ]; then
    echo "삭제 중: /etc/apt/sources.list.d/nginx.list"
    sudo rm -f /etc/apt/sources.list.d/nginx.list
  fi

  # NGINX 서명 키 파일 삭제
  if [ -f /usr/share/keyrings/nginx-archive-keyring.gpg ]; then
    echo "삭제 중: /usr/share/keyrings/nginx-archive-keyring.gpg"
    sudo rm -f /usr/share/keyrings/nginx-archive-keyring.gpg
  fi
  
  echo "패키지 목록을 다시 업데이트합니다..."
  sudo apt-get update >/dev/null 2>&1

  echo "롤백이 완료되었습니다. 시스템이 설치 시도 이전 상태로 복원되었습니다."
}

# 스크립트가 어떤 이유로든(오류 포함) 종료될 때 rollback 함수를 실행하도록 trap 설정
trap rollback EXIT

echo "Ubuntu 시스템에 NGINX 공식 패키지 설치를 시작합니다."
echo "========================================================"

# 이미 NGINX가 설치되어 있는지 확인
if command -v nginx &> /dev/null; then
    echo "NGINX가 이미 설치되어 있습니다. 설치를 중단합니다."
    # 정상 종료이므로 롤백이 실행되지 않도록 trap 해제
    trap - EXIT
    exit 0
fi

# 1단계: 필수 패키지 설치
echo "[1/6] 필수 패키지를 설치합니다..."
sudo apt-get update
sudo apt-get install -y curl gnupg2 ca-certificates lsb-release ubuntu-keyring
echo "필수 패키지 설치 완료."
echo ""

# 2단계: NGINX 공식 서명 키 가져오기 및 저장
echo "[2/6] NGINX 공식 서명 키를 가져옵니다..."
curl -fsSL https://nginx.org/keys/nginx_signing.key | gpg --dearmor \
    | sudo tee /usr/share/keyrings/nginx-archive-keyring.gpg >/dev/null
echo "서명 키를 저장했습니다."
echo ""

# 3단계: 서명 키 확인
echo "[3/6] 서명 키가 올바른지 확인합니다..."
FINGERPRINT_TO_CHECK="573BFD6B3D8FBC641079A6ABABF5BD827BD9BF62"
DOWNLOADED_FINGERPRINT=$(gpg --dry-run --quiet --no-keyring --import --import-options import-show /usr/share/keyrings/nginx-archive-keyring.gpg | grep "$FINGERPRINT_TO_CHECK")

if [ -z "$DOWNLOADED_FINGERPRINT" ]; then
    echo "오류: 다운로드한 NGINX 서명 키의 지문(fingerprint)이 공식 문서와 일치하지 않습니다."
    echo "설치 절차를 중단합니다."
    exit 1
else
    echo "키 지문(fingerprint)이 성공적으로 확인되었습니다."
fi
echo ""

# 4단계: NGINX APT 저장소 설정 (stable 버전)
echo "[4/6] NGINX stable 버전의 APT 저장소를 설정합니다..."
echo "deb [signed-by=/usr/share/keyrings/nginx-archive-keyring.gpg] \
http://nginx.org/packages/ubuntu `lsb_release -cs` nginx" \
    | sudo tee /etc/apt/sources.list.d/nginx.list
echo "APT 저장소 설정 완료."
echo ""

# 5단계: 저장소 Pinning 설정 (배포판 패키지보다 NGINX 저장소 패키지 우선)
echo "[5/6] 저장소 pinning을 설정하여 NGINX 공식 패키지를 우선합니다..."
echo -e "Package: *\nPin: origin nginx.org\nPin: release o=nginx\nPin-Priority: 900\n" \
    | sudo tee /etc/apt/preferences.d/99nginx
echo "Pinning 설정 완료."
echo ""

# 6단계: NGINX 설치
echo "[6/6] NGINX를 설치합니다..."
sudo apt-get update
sudo apt-get install -y nginx
echo "NGINX 설치 완료."
echo ""

echo "========================================================"
echo "✅ NGINX 설치가 성공적으로 완료되었습니다!"
echo ""
echo "아래 명령어를 사용하여 설치된 버전과 서비스 상태를 확인할 수 있습니다:"
echo "   nginx -v"
echo "   systemctl status nginx"
echo "========================================================"

# 모든 과정이 성공적으로 완료되었으므로, 스크립트가 정상 종료될 때 롤백 함수가 실행되지 않도록 trap을 해제합니다.
trap - EXIT

