#!/bin/sh

echo "Running pre-build script"

sed '/^export LIBZ_SYS_STATIC=1' -i /build.sh