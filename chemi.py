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

        e_c = 1.60217662*pow(10, -19)
        e_m = 9.10938356*pow(10, -31)
        planck = 6.62607004*pow(10, -34)
        Z = self.atomic_number
        ryd = 13.6057
        Z = Z - x
        E = -ryd*(math.pow(Z, 2)/math.pow(n, 2))

        return E

    def highest_energy(self, q=1):
        c = 0
        w = 0
        latest = 0
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

                        shell = (str(n+1) + s + str(m-l))
                        response['shells'][shell] = {}
                        response['shells'][shell]['energy'] = round(energy, 3)
                        response['shells'][shell]['number'] = self.shells[n][l][m]
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
    def __init__(self, name, entropy , enthalpy , constituents):
        self.constituents = constituents
        self.name = name
        self.entropy = float(entropy)
        self.enthalpy = float(enthalpy)
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
    def components(self):
        q = []
        for i in self.constituents:
            q.append(i.name)
        return q
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
            self.products.append(Compound("", 0, 0, self.compound))

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
            lnK = -(((self.enthalpy_change())/8.31) * (1/self.temperature)) + (self.entropy_change() / 8.31)
            return lnK
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

with open('data.csv', 'rb') as csvfile:
    datar = csv.reader(csvfile, delimiter='@', quotechar='|')
    for row in datar:
        preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[3])


with open('data2.csv', 'rb') as csvfile:
    datar = csv.reader(csvfile, delimiter='@', quotechar='|')
    for row in datar:
        preset_compound_data[row[0]] = Predefined_Compound(row[0], row[1], row[3])

################################################################################

if __name__ == "__main__":
    s = Reaction(300)
<<<<<<< HEAD
    a = Compound("CaCO3", 0, 0, [])
   # b = Compound("Al", 0, 0, [])


    s.reactants.append(a)
   # s.reactants.append(b)
=======
    a = Compound("NaCl", 0, 0, [])
    b = Compound("F2", 0, 0, [])


    s.reactants.append(a)
    s.reactants.append(a)
    s.reactants.append(b)
>>>>>>> 16bef7c04adce03cce3af8db09084eb7beea0a43


    s.predict()

    prod = s.return_products()
    react = s.return_reactants()

    print(output(react) + " -> " + output(prod))
    print("Gibbs Energy Change: " + str(s.gibbs_change()/1000) + "kJmol-1")
    print("Enthalpy Change: "+str(s.enthalpy_change()) + "kJmol-1")
    print("Entropy Change: "+str(s.entropy_change()) + "Jmol-1")
#TODO
