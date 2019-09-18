#!/usr/bin/env python3
import re 
with open ("data/Rsat_nonredundant_vertebrates_2017.tf", 'r') as matrices, open("data/matrixnames.txt", 'r') as names: 
	namelist = [name.rstrip("\n") for name in names.readlines()]
	mtxt = matrices.read()
	find = [re.compile('AC  '+name+"\n[^/]*?XX\n//", re.DOTALL) for name in namelist]
	#find = [re.compile("^AC.*"+name+".*^XX$\n//$", re.MULTILINE.DOTALL) for name in namelist]	 +"\n.*?XX\n//"
	found = [rgxp.search(mtxt).group(0) for rgxp in find if (rgxp.search(mtxt))]
	for match in found:
		print(match)


