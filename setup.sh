#!/usr/bin/env bash

# Setup Ubuntu
apt-get update && apt-get upgrade -y && apt-get install -y wget git gcc make tar build-essential libreadline-dev autotools-dev automake default-jre libncurses5-dev libncursesw5-dev

# install clann
git clone https://github.com/ChrisCreevey/clann.git
cd clann
# manually apply patches to code (not merged to master yet)
chmod +x install-sh
PATCH_LINE="                        while(string[i] != '(' && string[i] != ')' && string[i] != ',' && string[i] != ';' && string[i] != ':')"
sed -i "3236s/.*/${$PATCH_LINE}/" treecompare2.c
./configure
#make
# need to install manually due to some errors
gcc -DHAVE_CONFIG_H -I. -g -O2 -MT treecompare2.o -MD -MP -MF .deps/treecompare2.Tpo -c -o treecompare2.o treecompare2.c
mv -f .deps/treecompare2.Tpo .deps/treecompare2.Po
gcc -g -O2 -Wall -funroll-loops -o clann treecompare2.o -lreadline -lncurses -lm
mv clann /bin/clann
cd ..

echo "[+] Installed clann"

# install raxml
git clone https://github.com/stamatak/standard-RAxML.git
cd standard-RAxML
make -f Makefile.AVX.PTHREADS.gcc
rm *.o
mv raxmlHPC-PTHREADS-AVX /bin/raxml
cd ..

echo "[+] Installed RAxML"

# install ninja
curl http://nimbletwist.com/software/ninja/files/ninja.tgz > ninja.tgz
tar -xvf ninja.tgz
mv ninja_*/ninja /bin/ninja
mv ninja_*/Ninja.jar /bin/Ninja.jar

echo "[+] Installed ninja"

# Get Miniconda and make it the main Python interpreter
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh
chmod +x /miniconda.sh
/miniconda.sh -b -p /miniconda
rm /miniconda.sh

# install packages
/miniconda/bin/conda install -c bioconda -c conda-forge -c etetoolkit -c anaconda --file requirements.txt

echo "[+] Installed miniconda with packages"