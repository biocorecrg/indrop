FROM rocker/r-ver:3.4
MAINTAINER Viktor Petukhov "viktor.s.petuhov@ya.ru"

RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  less \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libz-dev \
  libxml2-dev \
  nano \
  pandoc

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools", "Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown", "withr"))'

#RUN useradd -m user
#USER user

#WORKDIR "/home/user"

RUN mkdir /project
WORKDIR "/project"

RUN \
  git clone --depth 1 https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/3.4

RUN \
  #echo '.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths()))' > ~/.Rprofile &&
  R -e 'chooseCRANmirror(ind=52); install.packages("ks", dependencies=c("Depends", "Imports", "LinkingTo"))' && \
  R -e 'devtools::install_local("/project/dropEst/dropestr/", dependencies=T)' && \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("dplyr", "ggplot2", "rlang"))'
RUN \
  mkdir -p /project/dropEst/build && \
  cd /project/dropEst/build && \
  cmake ../ && \
  make

COPY dropestr/dropReport.Rsc /project/dropEst/build/

ENV PATH=/project/dropEst/build:$PATH \
  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH
  #R_PROFILE=~/.Rprofile

ENTRYPOINT ["/bin/bash"]