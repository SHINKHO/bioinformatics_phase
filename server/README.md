1. 구동 환경별 서버 설정
    - production
        : nginx <==[uWSGI(worker)]==> [django app{user_app},django app{admin_app}]
    - local dev
        : [django apps ...]

        """
        > runserver시 auto-build ,보안 및 부하 문제 발생 
        """

2. DBMS : SQLite
    시범사업 및 닫힌(closed) 서버
    - Bioinformatics 툴 특성상 File-Based-Database를 사용
    - file-based 가 아닌 대부분의 db는 사용자, 프로젝트-초대,참여-권한 관련임.
    - django 는 ORM을 지원함
    ==> SQLite3 로 진행, 추후 Database 변경시 driver 하단만 바꿀수 있도록 ORM 기반 작성

3. 