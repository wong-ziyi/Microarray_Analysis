#! /bin/bash
cd ~/
sudo mkdir softwarepackages
cd softwarepackages
sudo curl -O https://cloud.r-project.org/bin/macosx/R-3.5.1.pkg
sudo curl -O https://www.python.org/ftp/python/2.7.15/python-2.7.15-macosx10.9.pkg
sudo curl -O https://root.cern.ch/download/root_v5.34.36.macosx64-10.11-clang70.tar.gz
sudo curl -O https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.OSX_x86_64.tar.gz
sudo curl -O https://jaist.dl.sourceforge.net/project/bowtie-bio/bowtie/1.2.2/bowtie-1.2.2-macos-x86_64.zip
sudo curl -O https://download1.rstudio.org/RStudio-1.1.456.dmg

sudo installer -pkg R-3.5.1.pkg -target /
sudo installer -pkg python-2.7.15-macosx10.9.pkg -target /
sudo pip3 install --upgrade pip
sudo pip3 install HTSeq

sudo gunzip -c root_v5.34.36.macosx64-10.11-clang70.tar.gz | sudo tar xopf -
sudo gunzip -c tophat-2.1.1.OSX_x86_64.tar.gz | sudo tar xopf -
sudo unzip -a bowtie-1.2.2-macos-x86_64.zip

sudo mv tophat-2.1.1.OSX_x86_64 tophat
sudo mv bowtie-1.2.2-macos-x86_64 bowtie

sudo rm -rf R-3.5.1.pkg
sudo rm -rf root_v5.34.36.macosx64-10.11-clang70.tar.gz
sudo rm -rf tophat-2.1.1.OSX_x86_64.tar.gz
sudo rm -rf bowtie-1.2.2-macos-x86_64.zip
sudo rm -rf python-2.7.15-macosx10.9.pkg

export PATH=$PATH:~/softwarepackages/tophat:~/softwarepackages/bowtie:~/softwarepackages/root/bin

sudo hdiutil attach RStudio-1.1.456.dmg
cd /Volumes/RStudio-1.1.456
sudo cp -rf Rstudio.app /Applications
sudo hdiutil detach /Volumes/RStudio-1.1.456
cd ~/
cd softwarepackages
sudo rm -rf RStudio-1.1.456.dmg
