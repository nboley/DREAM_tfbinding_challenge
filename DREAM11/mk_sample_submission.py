import os, sys
import gzip
import random

with gzip.open(sys.argv[1]) as fp:
    for line in fp:
        print line.strip() + "\t%.2e" % (random.random()**100)
