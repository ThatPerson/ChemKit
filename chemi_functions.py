import json
import re

with open('pt.json') as json_data:
	d = json.load(json_data)
	print (d['table'])
	json_data.close()
    
def get_chemical_object(c):
	for l in d['table']:
		
		if (l['small'] == c):
			return l
	print ("Chemical "+c+" not found!")
	return 0
    
def get_max_en(c):
	highest_val = 0
	highest_n = {}
	for i in c:
		if i['electronegativity'] > highest_val:
			highest_n = i
			highest_val = i['electronegativity']
	return highest_n
	
def get_min_en(c):
	lowest_val = 1000
	lowest_n = {}
	for i in c:
		if i['electronegativity'] < lowest_val and i['electronegativity'] != -1:
			lowest_n = i
			lowest_val = i['electronegativity']
	return lowest_n
    
def get_valency(c):
	outer_shell = c['electrons'][len(c['electrons']) - 1]
	valency = 1
	if (outer_shell < 8):
		if (outer_shell < 4):
			valency = outer_shell
		else:
			valency = 8 - outer_shell
	print(valency)   
	return valency 

chemicals_in_system = []

lo = input("> ")
comp = lo.split("+")
for i in comp: #This nets us a list of all the chemicals.
	print(i)
	chem = re.findall('[A-Z][a-z0-9]{0,3}', i)
	for p in chem:
		l = re.findall(r'\d+', p)
		if (len(l) == 1):
			num = int(l[0])
		else:
			num = 1
		q = re.sub(r'\d+', '', p)
		print (q)
		z = (get_chemical_object(q))
		if (z == 0):
			break
		for l in range(0, num):
			chemicals_in_system.append(z)
		
# The algorithm works on the following principle.
# The most reactive chemicals will react with the most reactive on the other end of the spectrum.
# We compute this with electronegativities - so if I have NaCl and HF then the highest EN (F) is paired with the lowest EN (Na). This gets NaF, and these two are removed from the chemicals_in_system. The process is then repeated.
chem_backup = chemicals_in_system
resultant_chemicals = []
while (len(chemicals_in_system) > 0):
	# First we get the highest electronegativity
	i = get_min_en(chemicals_in_system)
	print(i)
	chemicals_in_system.remove(i)
	# We assume maximum valency - so if oxygen is binding to carbon it will _always_ be a double bond (assuming carbon can fit it). 
	valency = get_valency(i)
	compound = [i]
	while valency > 0 and len(chemicals_in_system) > 0:
		q = get_max_en(chemicals_in_system)
		chemicals_in_system.remove(q)
		v = get_valency(q)
		v_left = 0
		if (valency < v):
			v_left = v-valency
			v = valency # Prevent excess shells being filled
		compound.append(q)
		valency = valency - v
		while (v_left > 0) and len(chemicals_in_system) > 0:
			y = get_min_en(chemicals_in_system) # Should probably do this stuff recursively, but I can't think straight.
			chemicals_in_system.remove(y)
			l = get_valency(y)
			l_left = 0
			if (v_left < l):
				# Then we cry.
				l = v_left
			compound.append(y)
			v_left = v_left - l
	resultant_chemicals.append(compound)
		


print(json.dumps(resultant_chemicals))
chemicals = []
output = lo + " -> "
first = 1
for i in resultant_chemicals:
	if (first == 0):
		output = output + " + "
	first = 0
	chemical = ""
	for p in i:
		chemical = chemical + p['small']
	print( chemical)
	chemicals.append(chemical)
	output = output + chemical
print (output)

#print (chemicals_in_system)


