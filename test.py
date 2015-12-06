import chemi

q = chemi.periodic_table["Na"].out(1)
print q["name"]
print q['shells']['3s0']
print q
s = chemi.Reaction(300)
q = chemi.Reaction(300)
s.reactants.append(chemi.Compound("LiBr", 0, 0, []))
s.reactants.append(chemi.Compound("LiBr", 0, 0, []))
s.reactants.append(chemi.Compound("Cl2", 0, 0, []))
s.predict()
print(chemi.output(s.return_reactants()) + " -> " + chemi.output(s.return_products()))
q.reactants = s.products
q.reactants[2] = chemi.Compound("F2", 0, 0, [])
q.predict()
print(chemi.output(q.return_reactants()) + " -> " + chemi.output(q.return_products()))
