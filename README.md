# ChemKit

ChemKit is a module I wrote quite a few years ago to help me write other chemistry related programs. I've written some documentation here.

```
import chemi
```

In the same directory. This gives access to the following classes;

### Element

This is a class which attempts to represent an element. In the current program, periodic_table is a dict of all elements so periodic_table['Na'] is an Element for Sodium.

Inside Element there are a number of variables;

| Name           | Contains |
|----------------|------------------------|
| name           | The long form name of an element (eg Sodium) |
| shells         | An array of arrays containing the predicted (or given) electronic shell configuration. Prediction is based on filling the lowest n, l first and filling m according to Hund's rules. So 1s, 2s, 2p, 3s, 3p, 3d, 4s - so it does not account for shielding leading to eg 4s being lower in energy than 3d. |
| small          | The symbol (eg Na) |
| position       | The group in the periodic table from 0 to 17 |
| molar          | The element's molar mass |
| atomic_number  | The element's atomic number |
| electronegativity | The element's Pauling electronegativity |
| electrons      | Arrangement of electrons by quantum number n |
| atomic_radius  | Radius in picometers (likely wrong as depends on definition of radius) |


There are also a number of functions

| Name           | Function |
|----------------|----------|
| get_shells()   | Fills shells from lowest n, l and m by Hund's rules |
| legacy_shell_energy(n,l,m,x) | Predicts energy, where x is the shielding, and ignoring any degeneracy from l & m - assumes Hydrogen like atom |
| shell_energy(n, l, m) | Wrong name, attempts to calculate effective nuclear charge using a bastardization of Slater's rules |
| highest_energy() | Uses shells to work out which is supposed to be the highest energy |
| outer_shell_energy() | Almost identical to highest_energy(), only this uses shell_energy() to give Zeff for the outer shell. So not energy. |
| out(dp) | Outputs JSON for an element giving properties, including predicted energy levels with energy and occupancy. dp refers to number of decimal places. |
| valency() | Returns predicted valency - from 0 or 8. Kind of useless on its own |

### Compound

Represents a compound made up of elements.

| Name         | Contains |
|--------------|----------|
| constituents | Array of constituent elements |
| name         | Chemical formula |
| entropy      | The entropy in J/Kmol |
| enthalpy     | Enthalpy in kJ/mol |
| molar        | Molar mass in atomic units |
| temp         | Temperature of model |

To assist this function, there are a number of functions to try and set up data. These are check_data(), get_new_data(), find_constituents(), find_name(), which try to obtain compound thermodynamic data from csv files (these do not contain values I've collected myself, however as I wrote this quite a few years ago I cannot remember where I got them from.), and to work out the constituents from name or name from constituents. Generally these run themselves.

More interesting functions are denoted below;

| Name        | Function |
|-------------|----------|
| melting_point() | Uses thermodynamic data from the attached csv files to identify the turning point for the reaction from the solid form to liquid form. Only works if the data is present, so fairly useless. |
| predict_mp_alg1() | Uses one of my very poor models to predict melting points, based on some very spurious assumptions for ionic compounds. Almost guaranteed to be wrong. |
| predict_mp_alg2() | See alg1. |
| predict_mp_alg3() | See alg1&2 |
| boiling_point()  | Uses thermodynamic data to identify boiling point - eg for H2O this gives 370.1K |
| sublimation_point() | Uses thermodynamic data to try and identify sublimation point |
| get_state(temperature) | Assuming all thermodynamic data is present this may be able to predict the state. But that's rare. |
| predict_bonding() | Tries to work out which elements would be bonded to which others. Best if you take it with a pinch of salt and assume it's talking about which bonds may be present. So for CH3OH it guesses O-H, C-O, C-H may be present. |
| mass() | Gives molar mass of compound |



Regarding the melting point predictions, I've run them for a few compounds along with the actual values below;

| Compound | Alg 1 | Alg 2 | Alg 3 | Actual (from Mathematica) |
|----------|-------|-------|-------|---------------------------|
| NaCl     | 1087.8| 1122.0| 0.26  | 1074.2                    |
| KCl      | 1097.4| |999.1| 0.00  | 1043.2                    |
| KF       | 1109.8| 1157.4| 0.31  | 1131.2                    |
| MgF2     | 965.1 | 865.0 | 0.20  | 987.2                     |
| H2O      | 1093.3| 3504.0| 0.26  | 273.2                     |

So generally for ionics Alg 2 tends to get the correct order, Alg 1 gets closest to the actual value, Alg 3 is just off doing its own thing, and for non ionics they all just collectively decide to not. They aren't based on any particular physical basis, it's just a bunch of spurious relationships.


### Reaction

Reaction tries to predict the outcome of a reaction between a number of compounds. It contains reactants and products (arrays of Compounds) and temperature. Running .predict() will populate products.

The assumption in this program is that very unelectronegative atoms will prefer to bond with very electronegative atoms. For example, Na will prefer to bond with F more than it would want to bond with Cl. And so the following occurs;

```
>>> s = chemi.Reaction(3)
>>> s.reactants.append(chemi.Compound("NaCl", 1, 0, 0, []))
>>> s.reactants.append(chemi.Compound("NaCl", 1, 0, 0, []))
>>> s.reactants.append(chemi.Compound("F2", 1, 0, 0, []))
>>> s.predict()
Data not found
[<chemi.Element instance at 0x7f0447f95dd0>, <chemi.Element instance at 0x7f0447f933b0>]
Data not found
Data not found
[<chemi.Element instance at 0x7f0447f95dd0>, <chemi.Element instance at 0x7f0447f933b0>]
Data not found
Data not found
[<chemi.Element instance at 0x7f044801acb0>, <chemi.Element instance at 0x7f044801acb0>]
Data not found
>>> s.products[0].name
u'NaF'
>>> s.products[1].name
u'NaF'
>>> s.products[2].name
u'Cl2'
>>>
```

So the reaction is predicted to be;

```
2NaCl + F2 -> 2NaF + Cl2
```

In this case, that is correct. This technique also works for eg

```
C3H8 + 5O2 -> 3CO2 + 4H2O
FeCl3 + Al -> AlCl3 + Fe
```

However it is very limited;

```
C2H6O + O -> C2H2 + 2H2O
```

Generally don't use it for organic molecules and don't trust it too much. It's less of an actual useful predictive tool and more of just an interesting thing I noticed a few years ago.

# Spectroscopy

In the Python directory there is a spectroscopy program. I wrote this a bit more recently, and it barely relates to ChemKit (other than using its periodic_table).

It uses a bunch of models to try to predict the frequencies of different rotational and vibrational features to produce a potential spectrum for a molecule. It gives output of a csv list, the first column of which is the frequency of the absorption (in cm-1) and the second is the intensity.
