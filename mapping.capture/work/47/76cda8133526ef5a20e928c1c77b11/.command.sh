#!/bin/bash -ue
printf 'Hello world!' | split -b 6 - chunk_
echo "ciao" > log/prova
