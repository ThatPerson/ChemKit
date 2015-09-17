import json
import re
import csv
import math

e_compounds = []
e_entropies = []
e_enthalpies = []

temp_k = 0

with open('data.csv', 'rb') as csvfile:
	datar = csv.reader(csvfile, delimiter=',', quotechar='|')
	for row in datar:
		e_compounds.append(row[0])
		e_entropies.append(row[3])
		e_enthalpies.append(row[1])

with open('data2.csv', 'rb') as csvfile:
	datar = csv.reader(csvfile, delimiter='@', quotechar='|')
	for row in datar:
		e_compounds.append(row[0])
		e_entropies.append(row[3])
		e_enthalpies.append(row[1])

preset_chemicals = [] # will be in the form ["varname", "chemical"]. Then I can just do the same thing as polyatomic.

compound = []
current = []
rounding = 5
#polyatomic = [["NO3", "Nr"]] # Polyatomic ions. Pretty shoddy way of dealing with it.
polyatomic = []
with open('pt.json') as json_data:
	d = json.load(json_data)

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
		if (i['electronegativity'] * ((i['position']/2)**2)) > highest_val:
			highest_n = i
			highest_val = (i['electronegativity'] * ((i['position']/2)**2))
	return highest_n
	
def get_min_en(c):
	lowest_val = 1000
	lowest_n = {}
	for i in c:
		if (i['electronegativity'] * ((i['position']/2)**2)) < lowest_val and (i['electronegativity'] * ((i['position']/2)**2)) != -1:
			lowest_n = i
			lowest_val = (i['electronegativity'] * ((i['position']/2)**2))
	return lowest_n
    
def get_valency(c):
	outer_shell = c['electrons'][len(c['electrons']) - 1]
	valency = 1
	if (outer_shell < 8):
		if (outer_shell < 4):
			valency = outer_shell
		else:
			valency = 8 - outer_shell

	return valency 

def chemical_to_name(c):
	chemicals = []
	number_of = []
	for i in c:
		if i['small'] in chemicals:
			number_of[chemicals.index(i['small'])] = number_of[chemicals.index(i['small'])] + 1
		else:
			chemicals.append(i['small'])
			number_of.append(1)
	chem = ""
	for i in range(0, len(chemicals)):
		chem = chem + chemicals[i]
		if (number_of[i] > 1):
			chem = chem + str(number_of[i])
			

	return chem
	
	

def get_next_set(mima, valency):
	global system_chemicals
	global compound
	global current
	#print ("Valency to fill is "+str(valency))
	while (valency > 0 and len(system_chemicals) > 0):
		if (mima == 1):
			plo = 0
			y = get_max_en(system_chemicals)
		else:
			plo = 1
			y = get_min_en(system_chemicals)

		if (y == {}):
			return 1
			
		system_chemicals.remove(y)
		v = get_valency(y)
		v_left = 0
		if (valency < v):
			v_left = v - valency
			v = valency
		compound.append(y)
		valency = valency - v
		
	       
		
		get_next_set(plo, v_left)
		
	return 1
    
    
	
# The algorithm works on the following principle.
# The most reactive chemicals will react with the most reactive on the other end of the spectrum.
# We compute this with electronegativities - so if I have NaCl and HF then the highest EN (F) is paired with the lowest EN (Na). This gets NaF, and these two are removed from the chemicals_in_system. The process is then repeated.
# Naturally, this algorithm doesn't necessarily get it right, but it does manage quite a few reactions well.
def get_resultant():
	global compound
	global system_chemicals
	chem_backup = system_chemicals
	resultant_chemicals = []		
       

	while (len(system_chemicals) > 0):
		# First we get the highest electronegativity
		i = get_min_en(system_chemicals)

		system_chemicals.remove(i)
		# We assume maximum valency - so if oxygen is binding to carbon it will _always_ be a double bond (assuming carbon can fit it). Could probably do it only as a single but it would produce a different outcome - I have no way to rank how good an outcome is. This would produce many different products.
					
		valency = get_valency(i)
		compound = [i]
		
		
		######## START HERE 
		get_next_set(1, valency)
		

				    
				    
		resultant_chemicals.append(compound)
	return resultant_chemicals
	
def find_chemical_system(c):
	c = c.replace(" ", "")
	global polyatomic
	global preset_chemicals
	for i in polyatomic:
	    c = c.replace(i[0], i[1])
	for p in preset_chemicals:
		c = c.replace(p[0], p[1])
	comp = c.split("+")
	chemicals_in_system = []
	for i in comp: #This nets us a list of all the chemicals.
		zx = i[0] 

		p = re.findall(r'\d+', zx)
		if (len(p) > 0):
			i = i[1:]
			lop = int(p[0])
			
		else:
			lop = 1
		chem = re.findall('[A-Z][a-z0-9]{0,3}', i)
		for p in chem:
			l = re.findall(r'\d+', p)
			if (len(l) == 1):
				num = int(l[0])
			else:
				num = 1
			q = re.sub(r'\d+', '', p)

			z = (get_chemical_object(q))
			if (z == 0):
				break
			for l in range(0, num):
				for ki in range(0, lop):
					chemicals_in_system.append(z)
	return chemicals_in_system
	
def get_mass(c):
	chem_sys = []
	chem_sys = find_chemical_system(c)
	mass = 0
	for i in chem_sys:
		mass = mass + i['molar']
	return mass

def get_compound_info(s):
	for i in range(0, len(e_compounds)):
		if (s == e_compounds[i]):
			return [float(e_entropies[i]), float(e_enthalpies[i])]
	return [0, 0]

def is_int(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False

chemicals_in_system = []
print ("Welcome to ChemKit (copyright 2015).")
print("This is primarily a reaction simulator.")
verbose = 0
lo = raw_input("> ")
while (lo != "exit"):
	system_chemicals = []
	qwo = lo.split(" ")
	
	zor = " ".join(qwo[1:])
	
	
	
	if (qwo[0] == "gibbs"):
		flip = 0
		reactants = []
		products = []
		reactants_entropy_total = 0
		reactants_enthalpy_total = 0
		products_entropy_total = 0
		products_enthalpy_total = 0
		t = 1
		co = ""
		for i in range(1, len(qwo)):
			if (qwo[i] == "->"):
				flip = 1
			else:
				if (qwo[i] != "+"):
					if (flip == 0):
						t = 1
						co = qwo[i]
						if (is_int(qwo[i][0])):
							t = int(qwo[i][0])
							co = qwo[i][1:]
						for p in range(0, t):
							reactants.append(co)
							x = (get_compound_info(co))
							reactants_entropy_total = reactants_entropy_total + x[0]
							reactants_enthalpy_total = reactants_enthalpy_total + (x[1]*1000)
					else:
						t = 1
						co = qwo[i]
						if (is_int(qwo[i][0])):
							t = int(qwo[i][0])
							co = qwo[i][1:]
						for p in range(0, t):
							products.append(co)
							x = (get_compound_info(co))
							products_entropy_total = products_entropy_total + x[0]
							products_enthalpy_total = products_enthalpy_total + (x[1]*1000)

		gibbss = (products_enthalpy_total - reactants_enthalpy_total) - (temp_k*(products_entropy_total - reactants_entropy_total))
		
		gibbs_products = products_enthalpy_total - (temp_k * products_entropy_total)
		gibbs_reactants = reactants_enthalpy_total - (temp_k * reactants_entropy_total)
		entropy_change = products_entropy_total - reactants_entropy_total
		enthalpy_change = products_enthalpy_total - reactants_enthalpy_total
		
		if (verbose == 1):
			print("Entropy Change of Reaction: "+str(entropy_change)+"Jmol-1K-1")
			print("Enthalpy Change of Reaction: "+str(enthalpy_change/1000)+"kJmol-1")
			print("Gibbs Free Energy at "+str(temp_k)+"K: "+str(gibbss/1000)+"kJmol-1")
		if (gibbss < 0):
			print("Will reaction go?: Yes")
		else:
			print("Will reaction go?: No")
		if (verbose == 1):
			print("Temperature: "+str(enthalpy_change/entropy_change)+"K")
			if (temp_k != 0):
				equilibrium = -(((enthalpy_change)/8.31) * (1/temp_k)) + (entropy_change / 8.31)
				print("ln K: "+str(equilibrium))
		


	
	
	if (qwo[0] == "resultant"):
		
		q = qwo[1:]
		
		chemicals_in_system = find_chemical_system(zor)
		system_chemicals = chemicals_in_system
		resultant_chemicals = get_resultant()

		chemicals = []
		number_of_c = []
		for ii in preset_chemicals:
			zor = zor.replace(ii[0], ii[1])

		output = zor + " -> "
		first = 1

		for i in resultant_chemicals:
			qwe = chemical_to_name(i)
			first = 0
			chemical = ""
	
			if (qwe in chemicals):
				number_of_c[chemicals.index(qwe)] = number_of_c[chemicals.index(qwe)] + 1
			else:
				chemicals.append(qwe)
				number_of_c.append(1)
		first = 1
	
		for i in range(0, len(chemicals)):
			if (first == 0):
				output = output + " + "
			first = 0
			if (number_of_c[i] > 1):
				output = output + str(number_of_c[i])
			output = output + chemicals[i]
		
		for i in polyatomic:
			output = output.replace(i[1], i[0])
		
		print (output)
		
		if (verbose == 1):
		    
		    
			# Eventually make it so that if verbose flag is set it prints out the formula in terms of the variables - so methane + 2oxygen -> etc.
			print ("reactants..")
			for i in q:
				if (i != '+'):
					l = i
					for poq in preset_chemicals:
						if (poq[0] == l):
							l = poq[1]
					print (l + ": " + (str(round(get_mass(i), rounding)) + "g/mol"))
			print ("products..")
			q_total = 0
			for i in range(0, len(chemicals)):
				q_total = q_total + (number_of_c[i] * get_mass(chemicals[i]))
			for i in range(0, len(chemicals)):
				qwewe = ''
				if (number_of_c[i] != 1):
					qwewe = str(number_of_c[i])
				print (qwewe + chemicals[i] + ": " + (str(round(number_of_c[i] * get_mass(chemicals[i]), rounding)) + "g/mol; ")+ str(number_of_c[i]*get_mass(chemicals[i])*100/q_total) + "%")
				
			
		
	if (qwo[0] == "mass"):
		print (str(round(get_mass(zor), 3)) + "g/mol")
	if (qwo[0] == "composition"):
		
		q = find_chemical_system(zor)
		total_mass = 0
		
		elements = []
		element_mass = []
		
		for i in q:
			total_mass = total_mass + i['molar']
			if (i['name'] in elements):
				element_mass[elements.index(i['name'])] = element_mass[elements.index(i['name'])] + i['molar']
			else:
				elements.append(i['name'])
				element_mass.append(i['molar'])
			
		if (total_mass > 0):
			for i in range(0, len(elements)):
				print(elements[i]+": "+str(round(100*(element_mass[i]/total_mass), rounding))+"%")
				
	if (qwo[0] == "set"):
		if (qwo[1] == "verbose"):
			verbose = 1
			print ("Verbose mode on")
		if (qwo[1] == "rounding"):
			if (len(qwo) > 2):
				rounding = int(qwo[2])
		if (qwo[1] == "temp"):
			temp_k = float(qwo[2])
	if (qwo[0] == "unset"):
		if (qwo[1] == "verbose"):
			verbose = 0
			print ("Verbose mode off")

	if (qwo[0] == "element"):
		if (len(qwo) > 1):
			for i in range(1, len(qwo)):
				p = get_chemical_object(qwo[i])
				if (p != 0):
					print ("=== "+p['name']+" ("+p['small']+") ===")
					print ("Atomic Number: "+str(p["number"]))
					print ("Atomic Mass: "+str(round(p["molar"], rounding)))
					print ("Electronegativity: "+str(round(p["electronegativity"], rounding)))
					
					#Perhaps put electron shells here?
					print ("")

	if (len(qwo) > 2):
		if (qwo[1] == "="):
			preset_chemicals.append([qwo[0], qwo[2]])
	if (qwo[0] == "help"):
		print "ChemSi is a chemical program - think Octave for chemicals. At the moment it is very basic."
		print "General commands are as follows. More will be added in the future;"
		print "set verbose - enables verbose mode (ie prints compositions in resultant)"
		print "mass [COMPOUND] - calculated the molecular mass of a compound."
		print "resultant [COMPOUND] + [COMPOUND] + ... - utilises an algorithm to predict the resultant chemical products of a reaction."
		print "composition [COMPOUND] - prints out the composition of a compound."
		print "gibbs [REACTANTS] -> [PRODUCTS] - prints out enthalpy change, entropy change, predicts optimum temperature and equilibrium constant"
	lo = raw_input("> ")
	
	
#print (chemicals_in_system)


