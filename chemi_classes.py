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

class Element:
    name = ""
    small = ""
    position = 0
    molar = 0
    atomic_number = 0
    electronegativity = 0
    lowest_energy_level = 0
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
    def __init__(self, name="", small="", position=0, molar=0, atomic_number=0, electronegativity=0):
        self.name = name
        self.small = small
        self.position = position
        self.molar = molar
        self.atomic_number = atomic_number
        self.electronegativity = electronegativity
        self.get_shells()
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

    def components(self):
        q = []
        for i in self.constituents:
            q.append(i.name)
        return q

    def find_consituents(self):
        name = self.name
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

# could have entropy and enthalpy values in list as they are now so they get added when the compound is created? also get rid of (g)/(l) etc before you run the split code - but maintain in the name (or a new phase variable?)

with open('pt.json') as json_data:
    d = json.load(json_data)
    for i in d['table']:
        periodic_table[i['small']] = Element(i['name'], i['small'], i['position'], i['molar'], i['number'], i['electronegativity'])
    json_data.close()

q = Compound("CaCO3")
q.find_consituents()
print q.components()
print q.find_name()
print q.composition()
