import sys
# Redirect print to a string buffer to capture the output for the final answer format
import io

old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Plan ---
# 1. Define monoisotopic masses for atoms and chemical moieties.
# 2. Calculate the initial mass of the unmodified glycan A2G2S2.
#    Composition: 4x GlcNAc, 5x Hexose (3 Man + 2 Gal), 2x Neu5Ac, linked together.
# 3. Calculate the total mass change from the chemical modifications:
#    a. Amidation: 2 sialic acids are amidated (COOH -> CONH2). This is equivalent to replacing -OH with -NH2.
#    b. Permethylation: All free -OH and -NH groups are methylated. We need to count the total number of sites.
# 4. Calculate the final mass of the singly sodiated ion [M+Na]+.
# 5. Print the step-by-step calculation and the final result.

# Step 1: Define monoisotopic masses
H = 1.007825
C = 12.000000
N = 14.003074
O = 15.994915
Na = 22.989770

# Mass of building blocks (as complete, non-residue molecules)
GLCNAC_MASS = (C*8) + (H*15) + (N*1) + (O*6) # 221.089938
HEXOSE_MASS = (C*6) + (H*12) + (O*6)         # 180.063388 (for Mannose and Galactose)
NEU5AC_MASS = (C*11) + (H*19) + (N*1) + (O*9) # 309.111121 (Sialic Acid)
H2O_MASS = (H*2) + O                         # 18.010565

# Step 2: Calculate initial mass of A2G2S2 glycan
num_glcnac = 4
num_hexose = 5 # 3 mannose + 2 galactose
num_neu5ac = 2
num_residues = num_glcnac + num_hexose + num_neu5ac
num_linkages = num_residues - 1

initial_mass = (num_glcnac * GLCNAC_MASS) + \
               (num_hexose * HEXOSE_MASS) + \
               (num_neu5ac * NEU5AC_MASS) - \
               (num_linkages * H2O_MASS)

# Step 3a: Calculate mass change from amidation
# Reaction is R-COOH + NH3 -> R-CONH2 + H2O.
# This is equivalent to replacing one -OH group with an -NH2 group.
mass_OH = O + H
mass_NH2 = N + (H * 2)
amidation_change_per_sia = mass_NH2 - mass_OH
total_amidation_change = 2 * amidation_change_per_sia # Two sialic acids

# Step 3b: Calculate mass change from permethylation
# Permethylation replaces H with CH3 on all free OH and NH groups.
# The mass added per site is Mass(CH3) - Mass(H) = Mass(CH2)
methylation_add_mass = C + (H*2) # 14.01565

# Counting methylation sites on the amidated glycan:
# Man3GlcNAc2 core:
#  - Reducing end GlcNAc: Me at N2, O1, O3, O6 -> 4 sites
#  - Core GlcNAc: Me at N2, O3, O6 -> 3 sites
#  - Man(b1-4)Man: Me at O2, O4 -> 2 sites (This is the branching Man)
#  - 2x Antennae Man: Me at C3, C4, C6 -> 3 sites/Man * 2 -> 6 sites
# Core total = 4 + 3 + 2 + 6 = 15 sites
# 2x Antennae GlcNAc: Me at N2, O3, O6 -> 3 sites/GlcNAc * 2 -> 6 sites
# 2x Antennae Gal: Me at O2, O4, and (O3 or O6) -> 3 sites/Gal * 2 -> 6 sites
# 2x Amidated Neu5Ac:
#  - Me at 4 -OH groups (C4, C7, C8, C9) -> 4 sites
#  - Me at N-acetyl N-H -> 1 site
#  - Me at both new amide N-H groups (-CONH2) -> 2 sites
#  - Total per Sia = 7 sites * 2 -> 14 sites
num_methylation_sites = 15 + 6 + 6 + 14
total_methylation_change = num_methylation_sites * methylation_add_mass

# Step 4: Calculate final m/z
final_mz = initial_mass + total_amidation_change + total_methylation_change + Na

# Step 5: Print the results clearly
print("All three glycans (A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2) will have the same mass after the described derivatization.")
print("\nThe expected m/z for the singly sodiated ion [M+Na]+ is calculated as follows:\n")
print("1. Mass of native A2G2S2 glycan:")
print(f"   ({num_glcnac}*m(GlcNAc) + {num_hexose}*m(Hex) + {num_neu5ac}*m(Neu5Ac)) - {num_linkages}*m(H2O) = {initial_mass:.4f} Da")
print("\n2. Mass change from derivatization:")
print(f"   Amidation of 2 sialic acids: 2 * ({amidation_change_per_sia:.4f}) = {total_amidation_change:.4f} Da")
print(f"   Permethylation at {num_methylation_sites} sites: {num_methylation_sites} * {methylation_add_mass:.4f} = {total_methylation_change:.4f} Da")
print("\n3. Addition of sodium ion for detection:")
print(f"   Mass of sodium adduct [M+Na]+: +{Na:.4f} Da")
print("\n--- Final Calculation ---")
print("m/z = Mass_native + Change_amidation + Change_methylation + Mass_Na")
print(f"m/z = {initial_mass:.4f} + ({total_amidation_change:.4f}) + {total_methylation_change:.4f} + {Na:.4f}")
print(f"m/z = {final_mz:.4f}")

# Restore stdout and get the captured output
sys.stdout = old_stdout
output_str = captured_output.getvalue()
print(output_str)

# Extract final number for the answer format
final_answer_val = f"{final_mz:.4f}"
print(f"<<<{final_answer_val}>>>")