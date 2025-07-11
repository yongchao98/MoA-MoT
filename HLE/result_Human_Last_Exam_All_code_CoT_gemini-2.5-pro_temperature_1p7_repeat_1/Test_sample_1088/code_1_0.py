# The 6-311G** basis set's primitive count depends on the atom.
# Let's calculate the total for a common molecule, methane (CH4).

# Define the number of atoms in methane
num_carbon_atoms = 1
num_hydrogen_atoms = 4

# --- Calculation for one Carbon (non-Hydrogen) atom ---
# Core primitives from the "6-" part
c_core_primitives = 6
# Valence primitives from the "-311" part (3+1+1)
c_valence_primitives = 3 + 1 + 1
# Polarization primitives from the first "*" (d-functions)
c_polarization_primitives = 1
# Total for one Carbon atom
c_total = c_core_primitives + c_valence_primitives + c_polarization_primitives

# --- Calculation for one Hydrogen atom ---
# Core primitives (Hydrogen has no core)
h_core_primitives = 0
# Valence primitives from the "-311" part (3+1+1)
h_valence_primitives = 3 + 1 + 1
# Polarization primitives from the second "*" (p-functions)
h_polarization_primitives = 1
# Total for one Hydrogen atom
h_total = h_core_primitives + h_valence_primitives + h_polarization_primitives

# --- Calculate the total for the methane (CH4) molecule ---
total_molecule_primitives = (num_carbon_atoms * c_total) + (num_hydrogen_atoms * h_total)

# --- Print the final breakdown and equation ---
print("Calculation for a Methane (CH4) molecule with the 6-311G** basis set:\n")
print(f"Primitives per Carbon atom = {c_core_primitives} (core) + {c_valence_primitives} (valence) + {c_polarization_primitives} (polarization) = {c_total}")
print(f"Primitives per Hydrogen atom = {h_core_primitives} (core) + {h_valence_primitives} (valence) + {h_polarization_primitives} (polarization) = {h_total}\n")
print("Total primitive Gaussians for CH4:")
print(f"({num_carbon_atoms} * {c_total}) + ({num_hydrogen_atoms} * {h_total}) = {total_molecule_primitives}")
