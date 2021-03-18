FROM coinstacteam/coinstac-base:centos7

WORKDIR /computation
#-------------------------------------------------
# Install Miniconda, and set up Python environment
#-------------------------------------------------
# Install required libraries
RUN yum install -y -q gcc libXext.x86_64 libXt.x86_64 \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/*
ENV PATH=/opt/miniconda/envs/default/bin:$PATH
RUN curl -ssL -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash miniconda.sh -b -p /opt/miniconda \
    && rm -f miniconda.sh \
    && /opt/miniconda/bin/conda update -n base conda \
    && /opt/miniconda/bin/conda config --add channels conda-forge \
    && /opt/miniconda/bin/conda create -y -q -n default python=3.7 \
    && rm -rf /opt/miniconda/[!envs]*

RUN npm install -g bids-validator

#Install other python packages
RUN pip install med2image
RUN pip install pillow tornado==5.0.2


#----------------------
# Install MCR and SPM12
#----------------------
# Install MATLAB Compiler Runtime
WORKDIR /opt
RUN curl -sSL -o mcr.zip https://ssd.mathworks.com/supportfiles/downloads/R2018b/deployment_files/R2018b/installers/glnxa64/MCR_R2018b_glnxa64_installer.zip \
    && unzip -q mcr.zip -d mcrtmp \
    && mcrtmp/install -destinationFolder /opt/mcr -mode silent -agreeToLicense yes \
    && rm -rf mcrtmp mcr.zip /tmp/*

 # Install standalone SPM
WORKDIR /opt
RUN curl -sSL -o spm.zip http://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/soon_gone/spm12_latest_Linux_R2018b.zip \
    && unzip -q spm.zip \
    && unzip /opt/spm12/spm12.ctf -d /opt/spm12/ \
    && rm -rf spm.zip
ENV MATLABCMD=/opt/mcr/v*/toolbox/matlab \
    SPMMCRCMD="/opt/spm*/run_spm*.sh /opt/mcr/v*/ script" \
    FORCE_SPMMCR=1 \
    LD_LIBRARY_PATH=/opt/mcr/v*/runtime/glnxa64:/opt/mcr/v*/bin/glnxa64:/opt/mcr/v*/sys/os/glnxa64:$LD_LIBRARY_PATH

RUN chmod -R 777 /opt/spm12

#Remove user warning from dicom init file
#RUN sed -i '53d' /opt/miniconda/envs/default/lib/python3.7/site-packages/dicom/__init__.py
#RUN sed -i '6d' /opt/miniconda/envs/default/lib/python3.7/site-packages/bids/grabbids/__init__.py

ADD . /computation
RUN pip install --no-cache-dir -r /computation/requirements.txt

# allow all inputs on vbm interface
COPY preprocess.py /opt/miniconda/envs/default/lib/python3.7/site-packages/nipype/interfaces/spm/
