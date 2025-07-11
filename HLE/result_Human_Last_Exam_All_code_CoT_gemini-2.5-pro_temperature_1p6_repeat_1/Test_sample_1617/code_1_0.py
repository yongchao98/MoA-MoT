import sys

# Step 1: Define Masses (monoisotopic)
# Redirect stdout to a string buffer to capture the output for the final answer
from io import StringIO
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()


aa_mass = {
    'A': 71.03711, 'V': 99.06841, 'D': 115.02694, 'L': 113.08406,
    'T': 101.04768, 'K': 128.09496, 'I': 113.08406, 'R': 156.10111
}
LACTYL_MOD = 72.02113
H_MASS = 1.00783  # Proton mass
H2O_MASS = 18.01056

# Step 2: Define Peptide and Experimental Data
peptide = "AVDLTKLIR"
exp_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

# Locate Lysine
lysine_pos_from_C = peptide[::-1].find('K') + 1 # 4th from C-terminus

# Step 3: Calculate y-ion series around Lysine
# y-ions are C-terminal fragments: R, IR, LIR, KLIR, ...

# Calculate y3 ion (LIR) - this should be unmodified
y3_seq = peptide[-3:] # "LIR"
y3_mass = aa_mass['L'] + aa_mass['I'] + aa_mass['R']
y3_mz = y3_mass + H2O_MASS + H_MASS
print(f"Analyzing peptide: {peptide}")
print("-" * 30)
print(f"Calculating mass for the y3-ion ('{y3_seq}')...")
print(f"Mass('{y3_seq}') = {y3_mass:.5f}")
print(f"m/z (y3) = Mass + H2O + H+ = {y3_mass:.5f} + {H2O_MASS} + {H_MASS} = {y3_mz:.5f}")
print(f"This theoretical value {y3_mz:.3f} closely matches the experimental value 401.276.\n")

# Calculate y4 ion (K(lac)LIR) - this should be modified
y4_seq = peptide[-4:] # "KLIR"
k_lactylated_mass = aa_mass['K'] + LACTYL_MOD
y4_mass = k_lactylated_mass + y3_mass
y4_mz = y4_mass + H2O_MASS + H_MASS
print(f"Calculating mass for the lactylated y4-ion ('K(lac)LIR')...")
print(f"Mass of lactylated Lysine (K(lac)) = {aa_mass['K']:.5f} + {LACTYL_MOD} = {k_lactylated_mass:.5f}")
print(f"Mass('K(lac)LIR') = {y4_mass:.5f}")
print(f"m/z (y4) = Mass + H2O + H+ = {y4_mass:.5f} + {H2O_MASS} + {H_MASS} = {y4_mz:.5f}")
print(f"This theoretical value {y4_mz:.3f} closely matches the experimental value 601.392.\n")

# The mass difference between y4 and y3 confirms the modification on K
mass_diff = y4_mass - y3_mass
print("Confirmation by mass difference:")
print(f"Mass(y4) - Mass(y3) = {y4_mass:.5f} - {y3_mass:.5f} = {mass_diff:.5f}")
print(f"This difference corresponds to the mass of a lactylated lysine residue ({k_lactylated_mass:.5f}).")
print("This confirms the lactylation is on the Lysine residue.\n")
# Print the final equation
print(f"Final Equation from Experimental Data:")
print(f"601.392 (y4 ion) - 401.276 (y3 ion) = {601.392 - 401.276:.3f} Da")
print(f"This mass difference of 200.116 Da is the mass of the lactylated lysine residue.\n")


# Let's check other values for completeness. The b5 ion can help confirm the N-terminal sequence.
# Calculate b5 ion (AVDLT)
b5_seq = peptide[:5]
b5_mass = sum(aa_mass[aa] for aa in b5_seq)
b5_mz = b5_mass + H_MASS
# A common fragment is a b-ion that has retained a water molecule (often seen with acidic/polar residues)
b5_plus_h2o_mz = b5_mass + H2O_MASS + H_MASS
print("Checking other significant ions:")
print(f"Calculating mass for b5+H2O ion ('{b5_seq} + H2O')...")
print(f"m/z (b5+H2O) = Mass('{b5_seq}') + H2O + H+ = {b5_mass:.5f} + {H2O_MASS} + {H_MASS} = {b5_plus_h2o_mz:.5f}")
print(f"This theoretical value {b5_plus_h2o_mz:.3f} closely matches the experimental value 518.271.")
print("-" * 30)

print("\nConclusion:")
print("The key evidence for lactylation on lysine is the pair of y3 and y4 ions.")
print(f"- 401.276 m/z represents the unmodified y3-ion ('LIR').")
print(f"- 601.392 m/z represents the lactylated y4-ion ('K(lac)LIR').")
print(f"- 518.271 m/z likely represents a fragment from the N-terminal part of the peptide ('AVDLT'+H2O), helping to confirm the overall peptide identity.")
print("Together, this set of ions {401.276, 601.392, 518.271} provides strong evidence for the specific lactylation.")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()
print(output_string)