import json
import re
import csv
import math

port = 5655
boltzmann = 1.23 * math.pow(10, -23)
avogadro = 6.02*math.pow(10, 23)
planck = 6.626*math.pow(10, -34)
periodic_table = {}
preset_compound_data = {}
class Element:
	def __init__(self, name="", small="", position=0, molar=0, atomic_number=0, electronegativity=0, electrons=[]):
		self.name = name
		self.shells = [
					[[0]],
					[[0], [0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]],
					[[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]]
				] # Over the 5f orbital it is inaccurate as it starts to have other shells - like 5 has another shell after f. For others it is not entirely accurate as it doesn't value it on energy - so 4s is not filled first, and there isn't the split as in copper and chromium.
		self.small = small
		self.position = position
		self.molar = molar
		self.atomic_number = atomic_number
		self.electronegativity = electronegativity
		self.get_shells()
		self.electrons = electrons
		self.atomic_radius = 0

	def get_shells(self):
		'''
			This is a very ugly algorithm. It works by assuming that the lowest
			orbital has the lowest energy - so, 1s is lower than 2s, 2s than 2p,
			etc. This breaks down at 3d/4s and later on - so Copper and Chromium
			both break the rule. Regardless for the majority of reactions it
			will work - if an element doesn't then it is possible to set the
			shells manually with Element.shells.
		'''
		a_n = self.atomic_number
		current_n = 0
		current_l = 0
		current_m = 0
		while (a_n > 0):
			q = 0

			if (current_n > 7):
				break
			if (self.shells[current_n][current_l][current_m+current_l] >= 2):
				if (current_m < current_l):
					current_m = current_m + 1
				else:
					if (current_l < current_n and current_l < 5):
						current_l = current_l + 1
						current_m = -current_l
					else:
						current_n = current_n + 1
						current_l = 0
						current_m = 0
			elif (self.shells[current_n][current_l][current_m+current_l] == 1):
				for i in range(current_m+current_l, current_l + current_l + 1):
					if (self.shells[current_n][current_l][i] == 0 and q == 0):
						self.shells[current_n][current_l][i] = self.shells[current_n][current_l][i] + 1
						a_n = a_n - 1
						q = 1
						break
				if (q == 0):
					self.shells[current_n][current_l][current_m+current_l] = self.shells[current_n][current_l][current_m+current_l] + 1
					a_n = a_n - 1
			else:
				self.shells[current_n][current_l][current_m+current_l] = self.shells[current_n][current_l][current_m+current_l] + 1
				a_n = a_n - 1

		return 1

	def shell_energy(self, n, l, m, x): # x n is the shell - ie 1, 2, 3. l is the subshell; spdf, m is the orbital; 1,2,3, x is the number of electrons before
		# Uses a formula derived from Rydberg when setting the new energy level to infinity - ie it has escaped.

		Z = self.atomic_number
		ryd = 13.6057
		Z = Z - x
		E = -ryd*(math.pow(Z, 2)/math.pow(n, 2))

	#This is the Rydberg formula - 1/lambda = RZ^2(1/n^2 - 1/n0^2) with n0 set to infinity - so it goes to RZ^2/n^2. As 1/lambda is proportional to Energy (E=hc/Lambda) this gives the energy. As ryd here is effectively hc/lambda for Z=1, n=1 it means no hc correction is needed.
		# ok no not actually the energy, but E=hf=hv/lam so proportional.
		return E

	def highest_energy(self, q=1):
		c = 0
		w = 0
		energy = 0
		for n in range(0, len(self.shells)):
			for l in range(0, len(self.shells[n])):
				w = 0
				for m in range(0, len(self.shells[n][l])):

					if(self.shells[n][l][m] != 0):
						energy = self.shell_energy(n+1, l, m-l, c)
						w = w + self.shells[n][l][m]
				c = c + w

		return energy*q


	def out(self, dp):
		response = {}
		response['name'] = self.name
		response['small'] = self.small
		response['an'] = self.atomic_number
		response['molar'] = round(self.molar, dp)
		response['electronegativity'] = round(self.electronegativity, dp)

		response['shells'] = {}

		c = 0
		w = 0

		for n in range(0, len(self.shells)):
			for l in range(0, len(self.shells[n])):
				w = 0
				for m in range(0, len(self.shells[n][l])):

					if(self.shells[n][l][m] != 0):
						energy = self.shell_energy(n+1, l, m-l, c)
						s = ""
						if l == 0:s = "s"
						elif l == 1:s = "p"
						elif l == 2:s = "d"
						elif l == 3: s = "f"

						shell = (str(n+1) + s + str(m-l))
						response['shells'][shell] = {}
						response['shells'][shell]['energy'] = round(energy, dp)
						response['shells'][shell]['number'] = self.shells[n][l][m]

						w = w + self.shells[n][l][m]
				c = c + w
		return response

	def valency(self):
		outer_shell = self.electrons[len(self.electrons) - 1]
		valency = 1
		if (outer_shell < 8): # Surely this limits it?
			if (outer_shell < 4):
				valency = outer_shell
			else:
				valency = 8 - outer_shell
		return valency

	def cat_or_an(self):
		outer_shell = self.electrons[len(self.electrons) - 1]
		if (outer_shell < 8):
			if (outer_shell < 4):
				return 1
			else:
				return 0

		return -1

class Compound:
	def __init__(self, name, molar, entropy , enthalpy , constituents, temp=0):
		self.constituents = constituents
		self.name = name
		self.full_name = name
		self.entropy = float(entropy)
		self.enthalpy = float(enthalpy)
		self.molar = molar
		self.temp = temp

		if (self.entropy == 0):
			self.check_data()
		if (name != ""):
			self.find_constituents()
		if (name == ""):
			self.find_name()
	def check_data(self):
		try:
			self.entropy = float(preset_compound_data[self.name].entropy)
			self.enthalpy = float(preset_compound_data[self.name].enthalpy)

		except:
			print "Data not found"

	def get_new_data(self):
		self.entropy = float(preset_compound_data[self.full_name].entropy)
		self.enthalpy = float(preset_compound_data[self.full_name].enthalpy)
		return 1

	def components(self):
		q = []
		for i in self.constituents:
			q.append(i.name)
		return q

	def melting_point(self):
		"""# how to get a. atomic spacing. Add up all atomic spacing in compound and find average?
		spacing = []
		for i in range(0, len(self.constituents)):
			for j in range(0, len(self.constituents)):
				if i != j:
					spacing.append(self.constituents[i].atomic_radius + self.constituents[j].atomic_radius)


		a = sum(spacing)/len(spacing) # a is in pm
		print(a)
		debye = 0
		# cl is Lindemann constant, approximating at 0.22
		Tm = (4*math.pow(math.pi, 2)*self.mass()*0.22*math.pow(a*math.pow(10, -12),2))/(boltzmann)
		#Tm = 2*math.pi*self.mass()*math.pow(0.22, 2)*math.pow(debye, 2)*boltzmann/math.pow(planck, 2)

		return Tm
		# http://phycomp.technion.ac.il/~phsorkin/thesis/node4.html
		# Should only be for crystals and similar but given that generally a molecule exists as a crystal when solid (unless it's organic but this doesn't really work with organics)
		# Tm = 4pi^2*m*cl*a^2/kb"""

		# Calculating from first principles is a pain. For some compounds I guess I could get (l) and (s) and work out the Gibbs change.
		name = self.name
		name = re.sub(r"\(.\)", "", name)
		l = Reaction(0)
		l.reactants.append(Compound(name+"(s)", 1, 0, 0, []))
		l.products.append(Compound(name+"(l)", 1, 0, 0, []))
		return l.turning_point()

		#########https://en.wikipedia.org/wiki/Born%E2%80%93Land%C3%A9_equation - Lattice Enthalpy might be correlated in some way

		###Name,		LE,			MP,			BP,			N
		#LiF	,	 1036,		1118,		1949,	 306		   , 1089.13
		#LiCl,	 853 ,		878.1,		1655,	 467.3		 , 1076.73
		#LiBr,	 807 ,		823.2,		1538,	 1071.9		, 1022.98
		#LiI	,	 757 ,		719.2,		1450,	 1558.7		, 973.373

	def predict_mp_alg1(self): # Only has any slight chance of being anywhere close if it is an ionic compound. Generally only if it is an ionic compound with two atoms - NaCl, NaBr, but you may have some luck with things like FeCl2
		Zo = 0 # Number of protons in cation - ie 11 for Na
		No = 0 # Row of periodic table of cation - ie 3 for Na
		Zt = 0 # Number of protons in anion - ie 17 for Cl
		Nt = 0 # Row of periodic table of anion - ie 3 for Cl
		for i in self.constituents:
			if (i.cat_or_an() == 1):
				Zo = Zo + i.atomic_number
				No = len(i.electrons)
			elif (i.cat_or_an() == 0):
				Zt = Zt + i.atomic_number
				Nt = len(i.electrons)
		print(Zo)
		print(No)
		print(Zt)
		print(Nt)
		if (No == 0 or Nt == 0):
			return -1
		R = 13.6
		N = abs(R*(math.pow(Zo/No, 2) - math.pow(Zt/Nt, 2)))
		melting_point = -0.101334 * N + 1109.81

		#m			   = 1.06499
		#b			   = 240.659
		#https://en.wikipedia.org/wiki/Born%E2%80%93Land%C3%A9_equation
		#Average with current function.

		

		return melting_point

	def predict_mp_alg2(self): # Uses Born-Lande equation to calculate Lattice Enthalpy. LE is correlated with melting point of test compounds with coeff 0.93.
		melting_point = 0

		return melting_point

	def boiling_point(self):
		name = self.name
		name = re.sub(r"\(.\)", "", name)
		l = Reaction(0)
		l.reactants.append(Compound(name+"(l)", 1, 0, 0, []))
		l.products.append(Compound(name+"(g)", 1, 0, 0, []))
		return l.turning_point()

	def sublimation_point(self):
		name = self.name
		name = re.sub(r"\(.\)", "", name)
		l = Reaction(0)
		l.reactants.append(Compound(name+"(s)", 1, 0, 0, []))
		l.products.append(Compound(name+"(g)", 1, 0, 0, []))
		return l.turning_point()

	def get_state(self, temp=-34): # this algorithm is only as good as the enthalpy and entropy data it's based on. I'm hoping to use equations to calculate the MP and BP without enthalpy or entropy. If I can it may be possible to do some simulatenous equations to get the enthalpy and entropy back from it!
		if (temp == -34):
			temp = self.temp
		mp = self.melting_point()
		bp = self.boiling_point()
		sp = self.sublimation_point()
		if (bp == -1 or mp == -1):
			if (sp == -1):
				return "n"
			else:
				if (temp >= sp):
					return "g"
				else:
					return "s"

		if (temp >= bp):
			return "g"
		elif (bp > temp >= mp):
			return "l"
		else:
			return "s"


		return 1

	def find_constituents(self):
		name = self.name
		if (len(name) == 0):
			return 0
		name = re.sub(r'\([a-z]+\)', '', name)
		p = re.findall(r'\d+', name[0])
		if (len(p) > 0):
			name = name[1:]
			lop = int(p[0])
		else:
			lop = 1
		chem = re.findall('[A-Z][a-z0-9]{0,3}', name)
		for p in chem:
			l = re.findall(r'\d+', p)
			if (len(l) == 1):
				num = int(l[0])
			else:
				num = 1
			q = re.sub(r'\d+', '', p)
			try:
				for l in range(0, num):
					for ki in range(0, lop):
						self.constituents.append(periodic_table[q])
			except:
				print("Element "+q+" not found.")

	def find_name(self):
		chemicals = []
		number_of = []
		for i in self.constituents:
			if (i in chemicals):
				number_of[chemicals.index(i)] = number_of[chemicals.index(i)] + 1
			else:
				chemicals.append(i)
				number_of.append(1)
		chem = ""
		for i in range(0, len(chemicals)):
			chem = chem + chemicals[i].small
			if (number_of[i] > 1):
				chem = chem + str(number_of[i])
		self.name = chem
		return chem

	def mass(self):
		for i in self.constituents:
			self.molar = self.molar + i.molar
		return self.molar

	def composition(self):
		response = {}
		self.mass()
		element_mass = {}
		for i in self.constituents:
			try:
				element_mass[i.small] = element_mass[i.small] + i.molar
			except:
				element_mass[i.small] = i.molar
		for key in element_mass:
			response[key] = (element_mass[key] / self.molar) * 100


		return response

	def gibbs_energy(self, temperature):
		return self.enthalpy - (temperature * self.entropy)


class Predefined_Compound:
	def __init__(self, name, enthalpy, entropy):
		self.name = name
		self.entropy = entropy
		self.enthalpy = enthalpy

class Reaction:
	def __init__(self, temp):
		self.temperature = temp
		self.reactants = []
		self.products = []
		self.system = []
		self.current = []

	def limiting_factor(self):
		lowest_v = Compound("", 0, 0, 0, [])
		lowest_c = 1000
		for i in self.reactants:
			if i.molar < lowest_c:
				lowest_c = i.molar
				lowest_v = i
		return [lowest_v.name, lowest_c]

	def find_en(self, type="HIGH"):
		if (type == "HIGH"):
			current = 0
		else:
			current = 10000
		curr_element = Element()
		for i in self.system:
			if (i.highest_energy(-1) != -1):
				if (type == "HIGH"):
					if (i.highest_energy(-1) * pow((i.position / 2), 2) > current):
						curr_element = i
						current = i.highest_energy(-1) * pow((i.position / 2), 2)
				else:
					if (i.highest_energy(-1) * pow((i.position / 2), 2) < current):
						curr_element = i
						current = i.highest_energy(-1) * pow((i.position / 2), 2)
		return curr_element

	def get_next_set(self, mima, valency, last):
		while (valency > 0 and len(self.system) > 0):
			if (mima == "HIGH"):
				plo = "LOW"
			else:
				plo = "HIGH"
			y = self.find_en(mima)
			if (len(y.name) == 0):
				return 1
			self.system.remove(y)
			v = y.valency()
			v_left = 0
			if (valency < v):
				v_left = v - valency
				v = valency
			self.compound.append(y)
			valency = valency - v

			self.get_next_set(plo, v_left, y)

	def predict(self):
		'''
			This algorithm works on the principle that a reaction will _always_
			tend to its lowest energy state, and that those atoms most willing
			to give up their electrons will have them taken by those most
			willing to take them. It is therefore not always right - for example
				C2H6O + O -> C2H4O + H2O (oxidation of ethanol)
			would not work - instead it does
				C2H6O + O -> C2H2 + 2H2O
			which is incorrect. For quite a few reactions, such as those below
			it does work pretty well. It does this by matching the highest
			outside orbital energy with the lowest - so, if Li, H and F were
			placed together then F by far has the lowest, and Li has the highest
			. ChemSi would predict that LiF would form as this would be more
			stable - because if you had HF then Li is more reactive than H
			(Li in H2O?) so the H would be displaced to give LiF.

				C3H8 + 5O2 -> 3CO2 + 4H2O
				FeCl3 + Al -> AlCl3 + Fe
				2NaCl + F2 -> 2NaF + Cl2

			and many others - play around with it.
			For some reactions where the enthalpies and entropies are known -
			see data.csv and data2.csv it will also try to predict gibbs energy.

		'''
		self.products = []
		for i in self.reactants:
			for p in i.constituents:
				self.system.append(p)
		while (len(self.system) > 0):
			i = self.find_en("LOW")
			self.system.remove(i)
			val = i.valency()
			self.compound = [i]

			self.get_next_set("HIGH", val, i)
			sdsd = self.limiting_factor()
			self.products.append(Compound("", sdsd[1], 0, 0, self.compound))

	def return_products(self):
		q = {}
		for i in self.products:
			try:
				q[i.name] = q[i.name] + 1
			except:
				q[i.name] = 1
		return q

	def return_reactants(self):
		q = {}
		for i in self.reactants:
			try:
				q[i.name] = q[i.name] + 1
			except:
				q[i.name] = 1
		return q

	def entropy_change(self):
		# find total entropy of reactants
		entropy_change = 0
		for i in self.reactants:
			entropy_change = entropy_change - i.entropy
		for i in self.products:
			entropy_change = entropy_change + i.entropy
		return entropy_change

	def enthalpy_change(self):
		enthalpy_change = 0
		for i in self.reactants:
			enthalpy_change = enthalpy_change - i.enthalpy
		for i in self.products:
			enthalpy_change = enthalpy_change + i.enthalpy
		return enthalpy_change

	def gibbs_change(self):
		return (1000*self.enthalpy_change()) - (self.temperature * self.entropy_change())

	def equilibrium_point(self):
		if (self.temperature != 0):


			# Arrhenius equation
			K = math.exp((-self.gibbs_change())/(8.3145 * self.temperature))
			return K
		else:
			return -1

	def turning_point(self, g_e = 0):
		# G = H - TS
		# TS = H - G
		# T = (H-G)/S
		return (1000*self.enthalpy_change() - g_e) / self.entropy_change()



def output(q):
	resp = ""
	x = 0
	for key in q:
		if (x != 0):
			resp = resp + " + "
		else:
			x = 1
		if (q[key] > 1):
			resp = resp + str(q[key]) + key
		else:
			resp = resp + key
	return resp

# could have entropy and enthalpy values in list as they are now so they get added when the compound is created? also get rid of (g)/(l) etc before you run the split code - but maintain in the name (or a new phase variable?)

######################LOAD DATA#################################################

with open('pt.json') as json_data:
	d = json.load(json_data)
	for i in d['table']:
		periodic_table[i['small']] = Element(i['name'], i['small'], i['position'], i['molar'], i['number'], i['electronegativity'], i['electrons'])
	json_data.close()

with open('radii.csv', 'rb') as csvfile:
	datar = csv.reader(csvfile, delimiter=',', quotechar='@')
	for row in datar:
		radii = 10
		for i in range(3, 6):
			if (row[i] != "na" and re.search('[a-zA-Z]',row[i]) != False):
				radii = row[i] # some elements only have empirical, some only have calculated. If they don't have either then it will default to 10.
				print(row[i])
				break

		try:
			periodic_table[row[1]].atomic_radius = float(radii)
		except:
			5 + 5
# atomic number,symbol,name,empirical t,Calculated,van der Waals,Covalent (single bond),Covalent (triple bond),Metallic (data is from https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page), while wikipedia might not be the most reliable this page is referenced and is nicely presented.)


with open('data.csv', 'rb') as csvfile:
	datar = csv.reader(csvfile, delimiter=',', quotechar='@')
	for row in datar:
		preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[3])


with open('data2.csv', 'rb') as csvfile:
	datar = csv.reader(csvfile, delimiter=',', quotechar='@')
	for row in datar:
		preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[3])

with open('data3.csv', 'rb') as csvfile:
	datar = csv.reader(csvfile, delimiter=',', quotechar='@')
	for row in datar:
		preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[2])
		#https://www.chem.wisc.edu/deptfiles/genchem/netorial/modules/thermodynamics/table.htm

################################################################################

if __name__ == "__main__":
	q = periodic_table["Na"].out(1)
   # print q["name"]

	s = Reaction(300)
	s.reactants.append(Compound("CuCl2", 3, 0, 0, []))
	s.reactants.append(Compound("2NaOH", 3, 0, 0, []))

#	s.reactants.append(Compound("Na", 3, 0, 0, []))
	#s.reactants.append(Compound("Cl2", 1, 0, 0, []))

#FeCl3 + Al -> AlCl3 + Fe

	s.predict()

	print(output(s.return_reactants()) + " -> " + output(s.return_products()))
	x = s.limiting_factor()
   # print(x[0] + " is limiting")

	"""compounds = ["NaCl", "KCl", "NaF", "KF", "NaBr", "KBr", "LiCl", "LiF", "LiBr", "RbCl" ,"RbF", "RbBr"]
	for i in compounds:
		l = Compound(i, 1, 0, 0, [])
		print(i + ", "+ str(l.melting_point()))
"""

	p = Reaction(373)
	p.reactants.append(Compound("H2O(l)", 1, 0, 0, []))
	p.products.append(Compound("H2O(g)", 1, 0, 0, []))
	print(p.turning_point())


	s = Compound("H2O", 1, 0, 0, [])
	print(s.get_state(0))
	print(s.get_state(10))
	print(s.get_state(100))
	print(s.get_state(1000))
	print(s.get_state(10000))
	print(s.get_state(100000))

	d = Compound("NaCl", 1, 0, 0, [])
	print(d.predict_mp_alg1())

#C3H8 + 5O2 -> 3CO2 + 4H2O
#FeCl3 + Al -> AlCl3 + Fe
#2NaCl + F2 -> 2NaF + Cl2


#TODO
