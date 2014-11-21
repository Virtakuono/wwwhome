#!/usr/bin/python

import sys
import os

url = sys.argv[1]

print(url)

os.system('wget '+sys.argv[1])
filename = url[(url.rfind('/')+1):]

printers = ['psy195','psy342','psy344','tupajumi','psy198']
printers = printers + ['psu352','psu351','psu131']
printers = printers + ['maari2','maari3','maari5','maari6']
printers = printers + ['psy339a','psy339b']


filebody = filename[0:-4]
os.system('pdf2ps ' + filename)
for i in printers:
	print('printing to '+ i)
	os.system('lpr -P' + i + ' -o sides=two-sided-long-edge ' + filebody + '.ps')
os.system('rm ' + filebody + '.ps')
os.system('rm ' + filename)


