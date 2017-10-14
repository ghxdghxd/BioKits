#!/usr/bin/env python
#-*- coding:utf-8 -*-
try:
	import psyco
	psyco.full()
except ImportError:
	pass

import os

name = os.popen("cut -f 1 ESCC.cnas.tsv")
name = set(name.readlines())
d={}
for n in name:
    d[n.strip()]=[]

print len(d)


with open("ESCC.cnas.tsv") as a:
    for line in a.readlines():
        word = line.strip().split()
        print word
        d[word[0]].append(word[1])

with open("ESCC.cnas","w") as f:
    for i in d:
        d[i] = list(set(d[i]))
        f.write(i+"\t"+"\t".join(d[i])+"\n")
