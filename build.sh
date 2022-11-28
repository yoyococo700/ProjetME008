echo $1
param=$1
echo $param

if [$param=build];
then
    gcc main.c
fi
