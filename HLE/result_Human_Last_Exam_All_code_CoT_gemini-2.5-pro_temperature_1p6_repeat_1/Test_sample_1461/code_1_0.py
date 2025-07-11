def calculate_helix_ring_size(monomer_lengths):
    """
    Calculates the H-bond ring size for an i->i+3 bond.
    The formula is: RingSize = 4 + Sum of backbone lengths of intervening residues.
    For an i->i+3 bond, the intervening residues are i+1 and i+2.
    """
    # Base atoms in the H-bond ring formula
    base_atoms = 4
    # Sum of backbone lengths of the two intervening monomers
    sum_of_lengths = sum(monomer_lengths)
    # Total ring size
    ring_size = base_atoms + sum_of_lengths
    return ring_size

# --- Define Monomer Properties ---

# Alanine is an alpha-amino acid
L_ala = 3
print(f"Backbone length of Alanine (L_Ala): {L_ala}")

# --- Scenario 1: Standard Epsilon Amino Acid ---
# Epsilon-amino acid has 7 backbone atoms
L_eps_standard = 7
print(f"\nScenario 1: Assuming a standard Epsilon Amino Acid")
print(f"Backbone length of Epsilon monomer (L_Eps): {L_eps_standard}")

# The intervening monomers for an i->i+3 bond are one Alanine and one Epsilon
intervening_monomers_scen1 = [L_ala, L_eps_standard]

# Calculate the resulting ring size
helix_size_scen1 = calculate_helix_ring_size(intervening_monomers_scen1)
print(f"The ring size is calculated as: 4 + L_Ala + L_Eps")
print(f"Ring Size = 4 + {L_ala} + {L_eps_standard} = {helix_size_scen1}")
print(f"This corresponds to a {helix_size_scen1}-helix.")

# --- Scenario 2: Strain-induced compact Epsilon Amino Acid ---
# A common effect of cyclic strain is to favor a more compact conformation.
# We will model this by using the length of a shorter (gamma) amino acid.
# Gamma-amino acid has 5 backbone atoms.
L_eps_compact_eff = 5
print(f"\nScenario 2: Assuming the 'cyclically-strained' Epsilon monomer behaves like a compact Gamma-like monomer")
print(f"Effective backbone length of strained monomer (L_eff): {L_eps_compact_eff}")

# The intervening monomers are one Alanine and one compact Epsilon-like monomer
intervening_monomers_scen2 = [L_ala, L_eps_compact_eff]

# Calculate the resulting ring size
helix_size_scen2 = calculate_helix_ring_size(intervening_monomers_scen2)
print(f"The ring size is calculated as: 4 + L_Ala + L_eff")
print(f"Ring Size = 4 + {L_ala} + {L_eps_compact_eff} = {helix_size_scen2}")
print(f"This corresponds to a {helix_size_scen2}-helix.")

# --- Conclusion ---
print(f"\nThe two most plausible helical patterns correspond to a {helix_size_scen2}-helix and a {helix_size_scen1}-helix.")
print(f"The resulting pattern is therefore {helix_size_scen2}/{helix_size_scen1}.")
