g++ -o TestGff -I../CppUtilClasses/ TestGff.cpp

if [[ $SAMTOOLPATH == "" ]]; then
	echo "\$SAMTOOLPATH not defined"
	exit
fi

if [[ $CPPUTILCLASSES == "" ]]; then
	echo "\$CPPUTILCLASSES not specified"
	exit
fi

g++ -o TestBamUtil -I$CPPUTILCLASSES -I$SAMTOOLINCLUDEPATH -L$SAMTOOLLIBPATH -lm -lbam -lz TestBamUtil.cpp $SAMTOOLLIBPATH/libbam.a