#!/bin/sh

echo "Running pre-build script"
sed '/^export CXX.*/a export PATH="/opt/osxcross/target/bin:$PATH"\nexport LIBZ_SYS_STATIC=1' -i /build.sh