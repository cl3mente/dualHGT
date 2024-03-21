# syntax=docker/dockerfile:1

# Comments are provided throughout this file to help you get started.
# If you need more help, visit the Dockerfile reference guide at
# https://docs.docker.com/go/dockerfile-reference/

# Want to help us make this template better? Share your feedback here: https://forms.gle/ybq9Krt8jtBL3iCk7
FROM ubuntu:20.04
# ARG PYTHON_VERSION=3.10.2
# FROM python:${PYTHON_VERSION}-slim as base

# Prevents Python from writing pyc files.
# ENV PYTHONDONTWRITEBYTECODE=1
ENV DEBIAN_FRONTEND=noninteractive

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering.
#ENV PYTHONUNBUFFERED=1

WORKDIR /app
RUN mkdir /app/input

# Download python, pip and required dependencies.
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    gffread \
    sudo \
    g++ \
    make

RUN apt-get install -y python3

RUN ln -s /usr/bin/python3 /usr/bin/python
RUN apt install -y python3-pip
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,source=./requirements.txt,target=./requirements.txt \
    pip install -r requirements.txt

# configure the MatPlotLib tmp directory
ENV MPLCONFIGDIR=/tmp/matplotlib

# Include OrthoFinder in the image
RUN wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder.tar.gz && \
    tar xzvf OrthoFinder.tar.gz && \
    rm OrthoFinder.tar.gz

# Include KaKs_Calculator in the image
RUN wget https://ngdc.cncb.ac.cn/biocode/tools/1/releases/3.0/file/download?filename=KaKs_Calculator3.0.zip -O KaKs_Calculator3.0.zip && \
    unzip KaKs_Calculator3.0.zip && \
    rm KaKs_Calculator3.0.zip && \ 
    cd KaKs_Calculator3.0/src && \
    make

# Include paraAT in the image
RUN wget https://download.cncb.ac.cn/bigd/tools/ParaAT2.0.tar.gz && \
    tar -xf ParaAT2.0.tar.gz && \
    rm ParaAT2.0.tar.gz

# Switch to the non-privileged user to run the application.
# Run the application.
# Create a non-privileged user that the app will run under.
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN apt-get install locate -y
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

USER appuser

# Copy the source code into the container.
COPY . .

ENV PATH=$PATH:/app/KaKs_Calculator3.0/src:/app/ParaAT2.0:/app/OrthoFinder



VOLUME /app/input

# CMD python HGT.py -i /data/bioinf2023/PlantPath2023/genomeANDgff -OFr /data/bioinf2023/PlantPath2023/genomeANDgff/results/prot/f7b812/Results_Feb23  -v -nt 50
# ENTRYPOINT ["python","./HGT.py"]