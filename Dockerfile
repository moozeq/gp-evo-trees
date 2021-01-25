FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

COPY requirements.txt requirements.txt
COPY pipe.py pipe.py

RUN apt-get update && apt-get install -y curl
RUN bash -c "$(curl -fsSL https://raw.githubusercontent.com/moozeq/GP_EvoTrees/pipeline/setup.sh)"

WORKDIR /evopipe
