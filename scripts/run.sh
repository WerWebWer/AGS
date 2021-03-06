#!/bin/bash

# shellcheck disable=SC2164
cd build
ctest --extra-verbose --repeat-until-fail 1 --timeout 10 --build-and-test
# shellcheck disable=SC2103
cd ..

FILES_SEQ="build/bin/Sequence"
for file in $FILES_SEQ; do
        echo "--------------------------------"
        echo $(basename $file)
        echo "--------------------------------"
        valgrind --error-exitcode=1 --leak-check=full --show-leak-kinds=all ./$file
done

FILES_STD="build/bin/Thread"
for file in $FILES_STD; do
        echo "--------------------------------"
        echo $(basename $file)
        echo "--------------------------------"
        valgrind --error-exitcode=1 --leak-check=full --show-leak-kinds=all ./$file
done
