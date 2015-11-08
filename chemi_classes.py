import json
import re
import csv
import math
import SocketServer
import sys
import json
from BaseHTTPServer import BaseHTTPRequestHandler
port = 5655

periodic_table = {}
preset_compound_data = {}
class Element:
    name = ""
    small = ""
    position = 0
    molar = 0
    atomic_number = 0
    electronegativity = 0
    lowest_energy_level = 0
    electrons = []
    shells = [
                [[0]],
                [[0], [0, 0, 0]],
                [[0], [0, 0, 0], [0, 0, 0, 0, 0]],
                [[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]],
                [[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [[0], [0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

            ]
    def __init__(self, name="", small="", position=0, molar=0, atomic_number=0, electronegativity=0, electrons=[]):
        self.name = name
        self.small = small
        self.position = position
        self.molar = molar
        self.atomic_number = atomic_number
        self.electronegativity = electronegativity
        self.get_shells()
        self.electrons = electrons
    def get_shells(self):
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
                    if (current_l < current_n):
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

    def shell_energy(self, n, l, m, x):
        e_c = 1.60217662*pow(10, -19)
        e_m = 9.10938356*pow(10, -31)
        planck = 6.62607004*pow(10, -34)
        Z = self.atomic_number
        ryd = 13.6057
        Z = Z - x
        E = -ryd*(math.pow(Z, 2)/math.pow(n, 2))

        return E
    def out(self, dp):
        response = ""
        response = response +  ("=== "+self.name+" ("+self.small+") ==="  + "NEWLINE")
        response = response +  ("Atomic Number: "+str(self.atomic_number) + "NEWLINE")
        response = response +  ("Atomic Mass: "+str(round(self.molar, dp)) + "NEWLINE")
        response = response +  ("Electronegativity: "+str(round(self.electronegativity, dp)) + "NEWLINE")

        c = 0
        w = 0
        latest = 0
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
                        response = response +  (str(n+1) + s + str(m-l)+" - "+str(round(energy, 3))+"eV: "+str(self.shells[n][l][m]) + "NEWLINE")
                        latest = energy
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


class Compound:
    constituents = []
    name = ""
    entropy = 0
    enthalpy = 0
    molar = 0
    def __init__(self, name="", entropy=0, enthalpy=0, constituents=[]):
        self.constituents = constituents
        self.name = name
        self.entropy = entropy
        self.enthalpy = enthalpy
        self.check_data()

        if (self.name == ""):
            self.find_name()
        if (self.constituents == []):
            self.find_consituents()

    def check_data(self):
        try:
            self.entropy = preset_compound_data[self.name].entropy
            self.enthalpy = preset_compound_data[self.name].enthalpy
        except:
            print "Data not found"

    def components(self):
        q = []
        for i in self.constituents:
            q.append(i.name)
        return q

    def find_consituents(self):
        name = re.sub(r'\([a-z]+\)', '', self.name)
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

class Predefined_Compound:
    enthalpy = 0
    entropy = 0
    name = ""
    def __init__(self, name, enthalpy, entropy):
        self.name = name
        self.entropy = entropy
        self.enthalpy = enthalpy

class Reaction:
    temperature = 0
    reactants = []
    products = []
    system = []
    compound = [] # current working compound.
    def __init__(self, temp):
        self.temperature = temp

    def find_en(self, type="HIGH"):
        if (type == "HIGH"):
            current = 0
        else:
            current = 100
        curr_element = Element()
        for i in self.system:
            if (i.electronegativity != -1):
                if (type == "HIGH"):
                    if (i.electronegativity * pow((i.position / 2), 2) > current):
                        print "CHANGE +"
                        curr_element = i
                        current = i.electronegativity * pow((i.position / 2), 2)
                else:
                    if (i.electronegativity * pow((i.position / 2), 2) < current):

                        print "CHANGE -"
                        curr_element = i
                        current = i.electronegativity * pow((i.position / 2), 2)
        return curr_element

    def get_next_set(self, mima, valency, last): # Basically just copied from old ChemSi code
        while (valency > 0 and len(self.system) > 0):
            if (mima == "HIGH"):
                plo = "LOW"
            else:
                plo = "HIGH"
            y = self.find_en(mima)
            if (y == Element()):
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

    def find_resultants(self):
        for i in self.reactants:
            for p in i.constituents:
                self.system.append(p)
        while (len(self.system) > 0):
            i = self.find_en("HIGH")
            self.system.remove(i)
            val = i.valency()
            self.compound = [i]

            self.get_next_set("LOW", val, i)
            self.products.append(Compound("", 0, 0, self.compound))
#def __init__(self, name="", entropy=0, enthalpy=0, constituents=[]):

# could have entropy and enthalpy values in list as they are now so they get added when the compound is created? also get rid of (g)/(l) etc before you run the split code - but maintain in the name (or a new phase variable?)

######################LOAD DATA#################################################

with open('pt.json') as json_data:
    d = json.load(json_data)
    for i in d['table']:
        periodic_table[i['small']] = Element(i['name'], i['small'], i['position'], i['molar'], i['number'], i['electronegativity'], i['electrons'])
    json_data.close()

with open('data.csv', 'rb') as csvfile:
    datar = csv.reader(csvfile, delimiter='@', quotechar='|')
    for row in datar:
        preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[3])


with open('data2.csv', 'rb') as csvfile:
    datar = csv.reader(csvfile, delimiter='@', quotechar='|')
    for row in datar:
        preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[3])

################################################################################

s = Reaction(300)
s.reactants.append(Compound("NaCl"))
s.reactants.append(Compound("NaCl"))
s.reactants.append(Compound("F2")) # Make it so Reaction has a parser.


q = 1
xz = ""
zx = ""
for i in s.reactants:

    if (q == 1):
        q = 0
    else:
        zx = zx + " + "

    zx = zx + i.name

q = 1

s.find_resultants()

for i in s.products:

    if (q == 1):
        q = 0
    else:
        xz = xz + " + "

    xz = xz + i.name

print zx + " -> " + xz
