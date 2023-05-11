FROM ubuntu:22.04

LABEL version="1.0"
LABEL maintainer="Ksenija Mladenovic<mk223336m@student.etf.bg.ac.rs>"
LABEL maintainer="Emilija Nikic<ne223077m@student.etf.bg.ac.rs>"
LABEL maintainer="Aleksa Bogdanovic<ba223333m@student.etf.bg.ac.rs>"

RUN apt-get update && apt-get install -y python3-pip && apt-get install -y git && apt-get install -y libgl1-mesa-glx && apt-get install -y libglib2.0-0

COPY requirements.txt .
# RUN pip install SpaGCN
# RUN pip install opencv-python
RUN pip install --no-cache-dir -r requirements.txt

WORKDIR /opt
RUN git clone https://github.com/Master-Computational-Genomics-SpaGCN/SpaGCN.git

WORKDIR /opt/SpaGCN

COPY Dockerfile /opt