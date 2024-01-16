#!/bin/sh

echo "Running pre-build script"

sed '/^export CXX.*/aexport LIBZ_SYS_STATIC=1' -i /build.sh