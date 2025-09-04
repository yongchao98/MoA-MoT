import math

# --- Constraints from Problem Statement ---
initial_mass = 7.20  # g
moles_H2O_exp = 3.60 / 18.015
moles_O2_exp = 0.80 / (15.999 * 2)
moles_N2_exp = 2.24 / 22.4

print("--- Experimental Data ---")
print(f"Moles H2O: {moles_H2O_exp:.4f}")
print(f"Moles O2:  {moles_O2_exp:.4f}")
print(f"Moles N2:  {moles_N2_exp:.4f}")
print("-" * 25)

# --- Hypothesis: Define Candidate Salts and Reactions ---
# Salt A: Ammonium Nitrite (NH4NO2)
M_A = 14.007 * 2 + 1.008 * 4 + 15.999 * 2  # g/mol
atoms_A = 2 + 4 + 2
# Reaction A: NH4NO2 -> 1*N2 + 2*H2O + 0*O2
react_A = {'N2': 1, 'H2O': 2, 'O2': 0}

# Salt B: Ammonium Nitrate (NH4NO3)
M_B = 14.007 * 2 + 1.008 * 4 + 15.999 * 3  # g/mol
atoms_B = 2 + 4 + 3
# Reaction B: NH4NO3 -> 1*N2 + 2*H2O + 0.5*O2
react_B = {'N2': 1, 'H2O': 2, 'O2': 0.5}

candidates = {
    'A': {'formula': 'NH4NO2', 'M': M_A, 'atoms': atoms_A, 'reaction': react_A},
    'B': {'formula': 'NH4NO3', 'M': M_B, 'atoms': atoms_B, 'reaction': react_B}
}

# --- Verification ---
# Constraint 1: Equimolar mixture must produce the observed moles of products.
# We can solve for 'n' (moles of each salt) from any of the products.
# Let's use O2, as it only comes from salt B.
# moles_O2_exp = n * react_B['O2']
n_from_O2 = moles_O2_exp / candidates['B']['reaction']['O2']
print(f"Calculated moles of each salt (n) from O2: {n_from_O2:.4f}")

# Now, check if this 'n' is consistent with the other products.
# Total moles = n * (moles from A) + n * (moles from B)
moles_H2O_calc = n_from_O2 * (candidates['A']['reaction']['H2O'] + candidates['B']['reaction']['H2O'])
moles_N2_calc = n_from_O2 * (candidates['A']['reaction']['N2'] + candidates['B']['reaction']['N2'])

# Constraint 2: Product mole consistency
pass_H2O = math.isclose(moles_H2O_calc, moles_H2O_exp, rel_tol=1e-3)
pass_N2 = math.isclose(moles_N2_calc, moles_N2_exp, rel_tol=1e-3)
print(f"Constraint Check (Product Moles):")
print(f"  H2O: Calculated={moles_H2O_calc:.4f}, Experimental={moles_H2O_exp:.4f} -> {'Pass' if pass_H2O else 'Fail'}")
print(f"  N2:  Calculated={moles_N2_calc:.4f}, Experimental={moles_N2_exp:.4f} -> {'Pass' if pass_N2 else 'Fail'}")

# Constraint 3: Total mass consistency
# Total mass = n * M_A + n * M_B
total_mass_calc = n_from_O2 * (candidates['A']['M'] + candidates['B']['M'])
pass_mass = math.isclose(total_mass_calc, initial_mass, rel_tol=1e-3)
print(f"Constraint Check (Total Mass):")
print(f"  Mass: Calculated={total_mass_calc:.2f} g, Initial={initial_mass:.2f} g -> {'Pass' if pass_mass else 'Fail'}")

# --- Final Calculation ---
if pass_H2O and pass_N2 and pass_mass:
    total_atoms = candidates['A']['atoms'] + candidates['B']['atoms']
    print("\nAll constraints passed. The candidates are correct.")
    print(f"Salt A: {candidates['A']['formula']} ({candidates['A']['atoms']} atoms)")
    print(f"Salt B: {candidates['B']['formula']} ({candidates['B']['atoms']} atoms)")
    print(f"Total number of atoms in A + B = {total_atoms}")
else:
    print("\nHypothesis failed.")
