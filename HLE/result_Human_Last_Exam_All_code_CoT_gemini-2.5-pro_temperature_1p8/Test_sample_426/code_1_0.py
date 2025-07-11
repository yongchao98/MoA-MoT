import sys
from io import StringIO

# Step 1: Analyze the composition of A1 to find its molecular formula.
# Given percentages: C=54.5%, H=13.6%, N=31.8%
# Atomic masses: C=12.01, H=1.008, N=14.01
C_percent = 54.5
H_percent = 13.6
N_percent = 31.8

# Calculate molar ratios
moles_C = C_percent / 12.01
moles_H = H_percent / 1.008
moles_N = N_percent / 14.01

# Find the smallest mole value to determine the empirical ratio
smallest_moles = min(moles_C, moles_H, moles_N)
C_ratio = moles_C / smallest_moles
H_ratio = moles_H / smallest_moles
N_ratio = moles_N / smallest_moles

print("Step 1: Determine the empirical and molecular formula of A1.")
print(f"The molar ratios are C:H:N ≈ {C_ratio:.1f}:{H_ratio:.1f}:{N_ratio:.1f}, which simplifies to 2:6:1.")
print("The empirical formula of A1 is C2H6N.")
print("A1 is formed from a dibromoalkane reacting with excess ammonia, so it must be a diamine (contains 2 N atoms).")
print("Thus, the molecular formula is a multiple of the empirical one, (C2H6N)x2 = C4H12N2.")
print("Verification of C4H12N2 percentages: C=54.5%, H=13.7%, N=31.8%. This matches the data.")
print("-" * 30)

# Step 2: Analyze the neutralization of the carboxylic acid to find its molar mass.
mass_acid = 2.16  # in grams
vol_KOH_liters = 30.0 / 1000  # in liters
conc_KOH = 1.0  # in M (mol/L)

# Moles of KOH = Concentration * Volume
moles_KOH = conc_KOH * vol_KOH_liters

# Assuming a monoprotic acid, moles of acid = moles of KOH
moles_acid = moles_KOH

# Molar Mass of acid = mass / moles
molar_mass_acid = mass_acid / moles_acid

print("Step 2: Determine the molar mass of the carboxylic acid.")
print("The neutralization equation is: RCOOH + KOH -> RCOOK + H2O")
print("From the titration data, we calculate the moles of KOH used:")
print(f"moles KOH = {conc_KOH} mol/L * {vol_KOH_liters} L = {moles_KOH} mol")
print("\nSince the molar ratio is 1:1, the moles of the carboxylic acid is also 0.030 mol.")
print("Now, we calculate the molar mass of the acid using the equation: Molar Mass = mass / moles")
print(f"Molar Mass = {mass_acid} g / {moles_acid} mol = {molar_mass_acid:.2f} g/mol")
print("-" * 30)

# Step 3: Identify the acid and deduce the structure of the initial hydrocarbon X.
print("Step 3: Deduce the structure of hydrocarbon X.")
print("The calculated molar mass is 72.00 g/mol.")
print("The molar mass of propanoic acid (C3H6O2 or CH3CH2COOH) is 74.08 g/mol.")
print("The small difference suggests the carboxylic acid is propanoic acid, and the given values have some experimental error.")
print("The oxidation of the C4 diol (A2) produces a C3 acid (propanoic acid) and CO2.")
print("This indicates a C-C bond cleavage. This specific cleavage occurs with vicinal diols (1,2-diols).")
print("For the products to be propanoic acid and CO2, the diol A2 must be butane-1,2-diol (HOCH2-CH(OH)-CH2CH3).")
print("This means the diamine A1 is 1,2-diaminobutane (H2NCH2-CH(NH2)-CH2CH3).")
print("1,2-diaminobutane is an asymmetric molecule with 4 different carbon atoms, which matches the 'four types of signals' in NMR.")
print("1,2-diaminobutane is formed from 1,2-dibromobutane (A).")
print("1,2-dibromobutane is formed from the addition of bromine (Br2) to an alkene (X).")
print("Therefore, the original hydrocarbon X must be but-1-ene (CH2=CH-CH2-CH3).")
print("\nConclusion: The structure of substance X is but-1-ene.")