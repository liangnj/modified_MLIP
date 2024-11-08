#!/bin/bash

wget -O example-1.tar.xz https://www.dropbox.com/s/jxen51rnq8km1a1/example-1.tar.xz?dl=0
wget -O example-2.tar.xz https://www.dropbox.com/s/ms8kkxdjvunavm7/example-2.tar.xz?dl=0
wget -O example-3.tar.xz https://www.dropbox.com/s/e6rq1n4mk4wpfze/example-3.tar.xz?dl=0

tar -xJvf example-1.tar.xz
tar -xJvf example-2.tar.xz
tar -xJvf example-3.tar.xz
