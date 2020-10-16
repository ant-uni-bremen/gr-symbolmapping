FROM ubuntu:focal

RUN DEBIAN_FRONTEND="noninteractive" apt update && apt -y install tzdata
RUN apt-get update \
    && apt-get install -y build-essential git cmake \
    libvolk2-dev liborc-0.4-dev \
    python3-dev python3-distutils python3-pybind11

ENV PYTHONPATH=/usr/local/lib/python3/dist-packages:$PYTHONPATH
WORKDIR /opt/dbuild
