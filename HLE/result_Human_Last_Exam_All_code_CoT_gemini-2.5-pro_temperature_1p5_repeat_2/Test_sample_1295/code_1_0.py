import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Molecular Information ---
smiles_representation = "C(=N)(NCCC)N=NC(=N)(NCCC)"

# --- Verification Calculations ---
# Atoms
num_C = 8
num_H = 18
num_N = 6
num_O = 0

# Valence Electrons
val_e_C = 4
val_e_H = 1
val_e_N = 5
total_valence_electrons = num_C * val_e_C + num_H * val_e_H + num_N * val_e_N

# Molecular Weight
mw_C = 12.011
mw_H = 1.008
mw_N = 14.007
total_molecular_weight = num_C * mw_C + num_H * mw_H + num_N * mw_N

# --- Output ---
print(f"SMILES Representation: {smiles_representation}")
print("\n--- Verification ---")

# Print Valence Electron Calculation
print(f"Total Valence Electrons: {total_valence_electrons}")
print(f"{total_valence_electrons} = {num_C} * {val_e_C} + {num_H} * {val_e_H} + {num_N} * {val_e_N}")

# Print Molecular Weight Calculation
print(f"\nApproximate Molecular Weight: {total_molecular_weight:.3f}")
print(f"{total_molecular_weight:.3f} = {num_C} * {mw_C} + {num_H} * {mw_H} + {num_N} * {mw_N}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output and the final answer in the required format
print(output)
print(f"<<<{smiles_representation}>>>")