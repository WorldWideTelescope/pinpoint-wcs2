#
#	Get revision number from mercurial
#

echo "#define VERSION \"0.9.0\"" > version.h
echo "#define REVISION "\"`hg id -n`\" >> version.h
