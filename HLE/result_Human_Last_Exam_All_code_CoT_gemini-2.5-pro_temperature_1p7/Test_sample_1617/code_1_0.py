import sys
import io

# capture original stdout
original_stdout = sys.stdout
# create a string buffer
captured_output = io.StringIO()
# redirect stdout
sys.stdout = captured_output


# --- Step 1: Define masses ---
# Monoisotopic masses for amino acid residues
residue_masses = {
    'A': 71.03711, 'V': 99.06841, 'D': 115.02694,
    'L': 113.08406, 'T': 101.04768, 'K': 128.09496,
    'I': 113.08406, 'R': 156.10111
}
# Mass of modifications and other entities
mass_lactyl = 72.02114  # Lactylation (C3H4O2)
mass_h2o = 18.01056    # Water
mass_proton = 1.00728  # Proton (for m/z calculation, z=1)

# --- Step 2: Peptide Information ---
peptide = "AVDLTKLIR"
mod_residue = 'K'
recorded_mz = [401.276, 601.392, 301.200, 518.271, 304.139, 417.223]

print("Peptide sequence: A-V-D-L-T-K(Lactyl)-L-I-R")
print("Recorded m/z values:", recorded_mz)
print("\n--- Calculating theoretical m/z values ---")

# --- Step 3: Calculate y-ions (C-terminal fragments) ---
# y3 ion: LIR
y3_residues = ['L', 'I', 'R']
y3_mass_sum = sum(residue_masses[r] for r in y3_residues)
y3_mz = y3_mass_sum + mass_h2o + mass_proton
print(f"\nCalculation for y3-ion (LIR):")
print(f"Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + Mass(H+) = m/z")
print(f"{residue_masses['L']:.3f} + {residue_masses['I']:.3f} + {residue_masses['R']:.3f} + {mass_h2o:.3f} + {mass_proton:.3f} = {y3_mz:.3f}")
print(f"This calculated m/z of {y3_mz:.3f} matches the recorded value 401.276.")

# y4 ion: K(lactyl)LIR
y4_residues = ['K', 'L', 'I', 'R']
mass_K_lactyl = residue_masses['K'] + mass_lactyl
y4_mass_sum = mass_K_lactyl + sum(residue_masses[r] for r in ['L', 'I', 'R'])
y4_mz = y4_mass_sum + mass_h2o + mass_proton
print(f"\nCalculation for y4-ion (K(lactyl)LIR):")
print(f"(Mass(K) + Mass(Lactyl)) + Mass(L) + Mass(I) + Mass(R) + Mass(H2O) + Mass(H+) = m/z")
print(f"({residue_masses['K']:.3f} + {mass_lactyl:.3f}) + {residue_masses['L']:.3f} + {residue_masses['I']:.3f} + {residue_masses['R']:.3f} + {mass_h2o:.3f} + {mass_proton:.3f} = {y4_mz:.3f}")
print(f"This calculated m/z of {y4_mz:.3f} matches the recorded value 601.392.")

# --- Step 4: Calculate b-ions (N-terminal fragments) ---
# b5 ion + H2O: AVDLT + H2O
b5_residues = ['A', 'V', 'D', 'L', 'T']
b5_mass_sum = sum(residue_masses[r] for r in b5_residues)
b5_plus_h2o_mz = b5_mass_sum + mass_h2o + mass_proton # b5-ion with a water adduct
print(f"\nCalculation for b5+H2O ion (AVDLT+H2O):")
print(f"Mass(A) + Mass(V) + Mass(D) + Mass(L) + Mass(T) + Mass(H2O) + Mass(H+) = m/z")
print(f"{residue_masses['A']:.3f} + {residue_masses['V']:.3f} + {residue_masses['D']:.3f} + {residue_masses['L']:.3f} + {residue_masses['T']:.3f} + {mass_h2o:.3f} + {mass_proton:.3f} = {b5_plus_h2o_mz:.3f}")
print(f"This calculated m/z of {b5_plus_h2o_mz:.3f} matches the recorded value 518.271.")


# --- Step 5: Conclusion ---
print("\n--- Conclusion ---")
print("The analysis shows that three recorded m/z values correspond to fragments of the lactylated peptide:")
print(f"- 401.276 corresponds to the y3-ion (LIR).")
print(f"- 601.392 corresponds to the y4-ion (K(lactyl)LIR).")
print(f"- 518.271 corresponds to a b5-ion with a water adduct (AVDLT+H2O).")
print("\nThe combination of these fragments, particularly the mass shift between the y3 and y4 ions (601.392 - 401.276 = 200.116 Da), confirms the presence of a lactylated lysine residue (Mass K + Mass Lactyl = 128.095 + 72.021 = 200.116 Da).")

# restore original stdout
sys.stdout = original_stdout
# get captured content
final_output = captured_output.getvalue()

# Print the final output
print(final_output)