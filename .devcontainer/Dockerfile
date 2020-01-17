FROM ubuntu:18.04
# 修改源地址
RUN mv /etc/apt/sources.list /etc/apt/sources_backup.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic main restricted " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-updates main restricted " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic universe " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-updates universe " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic multiverse " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-updates multiverse " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-backports main restricted universe multiverse " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-security main restricted " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-security universe " >> /etc/apt/sources.list && \
  echo "deb http://mirrors.ustc.edu.cn/ubuntu/ bionic-security multiverse " >> /etc/apt/sources.list && \
  echo "deb http://archive.canonical.com/ubuntu bionic partner " >> /etc/apt/sources.list && \
  # 更新并安装国际化语言包以及中文字体
  apt-get update && \
  apt-get install -y locales && \
  locale-gen zh_CN && \
  locale-gen zh_CN.utf8 && \
  apt-get install -y ttf-wqy-microhei ttf-wqy-zenhei xfonts-wqy
# 设置系统语言环境为中文UTF-8
ENV LANG zh_CN.UTF-8  
ENV LANGUAGE zh_CN.UTF-8
ENV LC_ALL zh_CN.UTF-8
# 安装开发环境库
RUN apt install -y gcc gdb build-essential && \
  # libraries
  apt install -y libeigen3-dev libdlib-dev && \
  # boost 
  apt install -y libboost-graph-dev && \
  # build tools
  apt install -y cmake && \
  apt install -y ninja-build
# google test dependencies
RUN apt install -y libpng-dev && \
  apt install -y libjpeg-dev && \
  apt install -y libblas-dev && \
  apt install -y liblapack-dev && \
  apt install -y libsqlite3-dev