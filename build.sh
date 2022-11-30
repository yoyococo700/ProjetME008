
param=$1

if [ $param == "build" ]; then
    gcc 'Projet Maths.c'
fi

if [ $param == "2" ]; then
    gcc 'Projet Maths.c' -lm
    ./a.out
fi