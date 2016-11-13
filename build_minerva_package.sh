#!/bin/bash
reset

pushd .

# Remove Old versions
rm minerva_1.4.5.tar.gz
R CMD REMOVE --library=../install minerva

# build new version
R CMD build --no-build-vignettes --no-manual minerva

# install in our directory
R CMD INSTALL --library=../install minerva_1.4.5.tar.gz

popd
