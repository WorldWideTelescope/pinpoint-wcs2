#! /bin/bash
#
# In order to build a Conda package, we need a tarball online that is generated
# with the Git submodule code included, which GitHub's default-built tarballs do
# not include. This script will generate such a tarball.
#
# This machinery is currently NOT plugged into any automation ... but it could be!

tarbase=pinpoint-wcs2

version="$1"

set -uo pipefail

if [ -z "$version" ] ; then
    echo >&2 "error: \$1 must be the current version number (checked against Git tag)"
    exit 1
fi

cur_tag=$(git describe --tags --exact-match HEAD)

if [ -z "$cur_tag" -o v"$version" != "$cur_tag" ] ; then
    echo >&2 "error: trying to build tarball for version $version but HEAD is not tag v$version"
    exit 1
fi

if [ -n "$(git status --porcelain)" ]; then 
    echo >&2 "error: Git working tree does not appear to be clean"
    exit 1
fi

set -e

git submodule update --init --recursive
git archive --format tar -o base.tar HEAD
git submodule --quiet foreach --recursive 'git archive --format tar --prefix=$displaypath/ -o submodule.tar HEAD'
top=$(pwd)
git submodule --quiet foreach --recursive "tar --concatenate --file=${top}/base.tar submodule.tar"

version_tar="${tarbase}-${version}.tar"
mv base.tar "$version_tar"
rm -f "${version_tar}.gz"
gzip -9 "${version_tar}"
sha256sum "${version_tar}.gz"

git submodule --quiet foreach --recursive 'rm -f submodule.tar'
