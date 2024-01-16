#! /bin/bash

echo -e "\033[31m Generate a key pair\033[37m"
./NTRUSign -g MySignature
ls

echo -e "\033[31mCreate a random file\033[37m"
dd if=/dev/random of=file.test bs=1M count=$(($RANDOM % 100))
ls

echo -e "\033[31mSign the file\033[37m"
./NTRUSign -is MySignature_priv.asc -s file.test
ls

echo -e "\033[31mVerify the file\033[37m"
./NTRUSign -ip MySignature_pub.asc -v file.test

mv file.test file.test_bck

echo "--------------------------------------------------------"
echo -e "\033[31mModify the file\033[37m"
dd if=/dev/random of=file.test bs=1M count=$(($RANDOM % 100))
ls

echo -e "\033[31mVerify the file\033[37m"
./NTRUSign -ip MySignature_pub.asc -v file.test

echo "--------------------------------------------------------"

mv file.test_bck file.test

echo -e "\033[31mGenerate a bad key pair\033[37m"
./NTRUSign -g BadSignature
ls

echo -e "\033[31mSign the file\033[37m"
./NTRUSign -is BadSignature_priv.asc -s file.test
ls

echo -e "\033[31mVerify the file\033[37m"
./NTRUSign -ip MySignature_pub.asc -v file.test
