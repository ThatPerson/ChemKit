import chemi
import math

class SpectroscopySystem:
	temperature = 273
	electronic_energy_levels = []
	bond_length = 0
	vibrational_constant = 0
	anharmonicity_constant = 0
	rotational_constant = 0
	atom_a = 0
	atom_b = 0
	def reduced_mass(self):
		return (self.atom_a.molar * self.atom_b.molar) / (self.atom_a.molar + self.atom_b.molar)

	def __init__(self, ab, bb, energy_levels, beta, k, *args, **kwargs):
		self.electronic_energy_levels = energy_levels
		self.atom_a = ab
		self.atom_b = bb
		bl = kwargs.get('bl', ab.atomic_radius + bb.atomic_radius)
		self.bond_length = bl * pow(10, -12)
		we = kwargs.get('we', 5.0341*pow(10,22) * math.sqrt(k / ((ab.molar + bb.molar) * 1.66 * pow(10, -27)))) # in reciprocal centimeters
		xe = kwargs.get('xe', (1.0546 * pow(10, -34) * pow(beta, 2)) / (2 * ((ab.molar + bb.molar) * 1.66 * pow(10, -27)) * we))
		self.vibrational_constant = we
		self.anharmonicity_constant = xe
		I = self.reduced_mass() * 1.66*pow(10, -27) * pow(self.bond_length, 2)
		self.rotational_constant = kwargs.get('rot', (pow(6.626 * pow(10, -34), 1))/(8*pow(math.pi, 2) * (3*pow(10, 10)) * I))

		self.temperature = kwargs.get('t', 273)

	def transition(self, ei, ef, vi, vf, ji, jf):
		electronic = self.electronic_energy_levels[ef] - self.electronic_energy_levels[ei]
		Evi = (vi + 0.5)*self.vibrational_constant - (pow((vi + 0.5), 2) * self.vibrational_constant * self.anharmonicity_constant)
		Evf = (vf + 0.5)*self.vibrational_constant - (pow((vf + 0.5), 2) * self.vibrational_constant * self.anharmonicity_constant)
		Evg = (0 + 0.5)*self.vibrational_constant - (pow((0 + 0.5), 2) * self.vibrational_constant * self.anharmonicity_constant)
		vibrational = Evf - Evi
		Eji = self.rotational_constant * ji * (ji + 1)
		Ejf = self.rotational_constant * jf * (jf + 1)
		rotational = Ejf - Eji
		total = electronic+vibrational+rotational
		rot_population = ((2 * ji) + 1) * math.exp(- (Eji / (0.695 * self.temperature))) # of initial
		elec_population = math.exp(- (self.electronic_energy_levels[ei] / (0.695 * self.temperature)))
		vib_population = math.exp(- (Evi - Evg) / (0.695 * self.temperature))
		pop = rot_population * elec_population * vib_population
		return [total, pop]


x = SpectroscopySystem(chemi.periodic_table["H"], chemi.periodic_table["H"], [0], 3.8955*pow(10, -7), 2.53769*pow(10,-65), bl = 74)

lines = []
populations = []

for v_initial in range(0, 10):
	for v_final in range(v_initial - 3, v_initial + 3):
		if (v_final > 0 and v_final != v_initial):
			for j_initial in range(0, 10):
				if (j_initial > 0):
					xc = x.transition(0, 0, v_initial, v_final, j_initial, j_initial - 1)
					lines.append(xc[0])
					populations.append(xc[1])
				xc = x.transition(0, 0, v_initial, v_final, j_initial, j_initial + 1)
				lines.append(xc[0])
				populations.append(xc[1])

for i in range(0, len(lines)):
	print(str(lines[i]) + ", " + str(populations[i]))
