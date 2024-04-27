#!/bin/bash

# fastp installation
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

# diamond installation
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
chmod +x diamond
rm diamond-linux64.tar.gz
