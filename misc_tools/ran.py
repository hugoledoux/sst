#! /usr/bin/python

import sys
import random
import csv


with open('100.xyz', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=' ')
    for i in range(100):
        x = random.uniform(0, 100)
        y = random.uniform(0, 100)
        writer.writerow([x,y,0])




