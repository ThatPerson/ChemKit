import json
import re

compound = []

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
        #print ("Valency to fill is "+str(valency))
        while (valency > 0 and len(system_chemicals) > 0):
                if (mima == 1):
                        plo = 0
                        y = get_max_en(system_chemicals)
                else:
                        plo = 1
                        y = get_min_en(system_chemicals)
                system_chemicals.remove(y)
                v = get_valency(y)
                v_left = 0
                if (valency < v):
                        v_left = v - valency
                        v = valency
                compound.append(y)
                valency = valency - v
                #print ("Added "+y['small']+" to mixture. Remaining valency is "+str(valency))
                get_next_set(plo, v_left)
                
        return 1
    
    
        
# The algorithm works on the following principle.
# The most reactive chemicals will react with the most reactive on the other end of the spectrum.
# We compute this with electronegativities - so if I have NaCl and HF then the highest EN (F) is paired with the lowest EN (Na). This gets NaF, and these two are removed from the chemicals_in_system. The process is then repeated.
def get_resultant():
        chem_backup = system_chemicals
        resultant_chemicals = []                
        global compound
        global system_chemicals

        while (len(system_chemicals) > 0):
                # First we get the highest electronegativity
                i = get_min_en(system_chemicals)

                system_chemicals.remove(i)
                # We assume maximum valency - so if oxygen is binding to carbon it will _always_ be a double bond (assuming carbon can fit it). 
                valency = get_valency(i)
                compound = [i]
                
                
                ######## START HERE 
                get_next_set(1, valency)
                
                
                ''' while valency > 0 and len(system_chemicals) > 0:
                        q = get_max_en(system_chemicals)
                        system_chemicals.remove(q)
                        v = get_valency(q)
                        v_left = 0
                        if (valency < v):
                                v_left = v-valency
                                v = valency # Prevent excess shells being filled
                        compound.append(q)
                        valency = valency - v
                        while (v_left > 0) and len(system_chemicals) > 0:
                                y = get_min_en(system_chemicals) # Should probably do this stuff recursively, but I can't think straight.
                                system_chemicals.remove(y)
                                l = get_valency(y)
                                l_left = 0
                                if (v_left < l):
                                        # Then we cry.
                                        vp_left = l - v_left
                                        l = v_left
                                compound.append(y)
                                v_left = v_left - l
                                
                                while (vp_left > 0) and len(system_chemicals) > 0:
                                    yp = get_max_en(system_chemicals) # Should probably do this stuff recursively, but I can't think straight.
                                    system_chemicals.remove(yp)
                                    lp = get_valency(yp)
                                    lp_left = 0
                                    if (vp_left < l):
                                            # Then we cry.
                                            l = v_left
                                    compound.append(yp)
                                    vp_left = vp_left - l
                                    
                        '''            
                                    
                                    
                                    
                resultant_chemicals.append(compound)
        return resultant_chemicals
        
def find_chemical_system(c):
        c = c.replace(" ", "")
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

chemicals_in_system = []
print ("Welcome to ChemKit (copyright 2015).")
verbose = 0
lo = input("> ")
while (lo != "exit"):
        system_chemicals = []
        qwo = lo.split(" ")
        
        zor = " ".join(qwo[1:])
        if (qwo[0] == "resultant"):
                
                q = qwo[1:]
                
                chemicals_in_system = find_chemical_system(zor)
                system_chemicals = chemicals_in_system
                resultant_chemicals = get_resultant()

                chemicals = []
                number_of_c = []
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

                
                print (output)
                
                if (verbose == 1):
                        print ("reactants..")
                        for i in q:
                                if (i != '+'):
                                        print (i + ": " + (str(round(get_mass(i), 3)) + "g/mol"))
                        print ("products..")
                        q_total = 0
                        for i in range(0, len(chemicals)):
                                q_total = q_total + (number_of_c[i] * get_mass(chemicals[i]))
                        for i in range(0, len(chemicals)):
                                qwewe = ''
                                if (number_of_c[i] != 1):
                                        qwewe = str(number_of_c[i])
                                print (qwewe + chemicals[i] + ": " + (str(round(number_of_c[i] * get_mass(chemicals[i]), 3)) + "g/mol; ")+ str(number_of_c[i]*get_mass(chemicals[i])*100/q_total) + "%")
                                
                        
                
        if (qwo[0] == "mass"):
                

                print (str(round(get_mass(zor), 3)) + "g/mol")
        if (qwo[0] == "set"):
                if (qwo[1] == "verbose"):
                        verbose = 1
                        print ("Verbose mode on")

        if (qwo[0] == "unset"):
                if (qwo[1] == "verbose"):
                        verbose = 0
                        print ("Verbose mode off")
        lo = input("> ")
        
        
#print (chemicals_in_system)


