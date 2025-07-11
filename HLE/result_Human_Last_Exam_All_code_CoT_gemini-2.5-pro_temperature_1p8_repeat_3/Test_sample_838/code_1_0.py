import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

print("This script details the multi-step synthesis from [(3S)-3-bromobutyl]benzene.\n")

print("--- Reaction 1: Formation of Product A ---")
print("Starting Material: [(3S)-3-bromobutyl]benzene")
print("Reagents: Potassium tert-butoxide (t-BuOK) in cyclohexane/diethyl ether.")
print("Reaction Type: E2 Elimination.")
print("Explanation: Potassium tert-butoxide is a strong, bulky base. It favors Hofmann elimination, abstracting a proton from the least sterically hindered carbon (the terminal methyl group). This forms the less substituted alkene.")
print("Product A is 4-phenylbut-1-ene. The original chiral center is destroyed, making product A achiral.")
print("-" * 20)

print("\n--- Reaction 2: Formation of Product B ---")
print("Starting Material: Product A (4-phenylbut-1-ene)")
print("Reagents: 1. Borane in THF (BH3·THF), 2. Hydrogen peroxide (H2O2) and sodium hydroxide (NaOH).")
print("Reaction Type: Hydroboration-Oxidation.")
print("Explanation: This is an anti-Markovnikov addition of water across the double bond. The hydroxyl group (-OH) is added to the terminal, less substituted carbon of the alkene.")
print("Product B is 4-phenylbutan-1-ol. This product is achiral.")
print("-" * 20)

print("\n--- Reaction 3: Formation of Final Product C ---")
print("Starting Material: Product B (4-phenylbutan-1-ol)")
print("Reagent: Phosphorous tribromide (PBr3).")
print("Reaction Type: Substitution (SN2).")
print("Explanation: PBr3 is a reagent that converts primary alcohols into primary alkyl bromides by replacing the -OH group with a -Br atom.")
print("The final product, C, is 1-bromo-4-phenylbutane.")
print("-" * 20)

print("\n--- Identity of Final Product C ---")
final_product_name = "1-bromo-4-phenylbutane"
print(f"IUPAC Name: {final_product_name}")

print("\nChirality Explanation:")
print("The original molecule, [(3S)-3-bromobutyl]benzene, was chiral. However, the chiral center at carbon-3 was destroyed during the first step (E2 elimination), which resulted in the formation of an achiral alkene. No new stereocenters were formed in the subsequent reactions. Therefore, the final product C, 1-bromo-4-phenylbutane, is achiral.")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)

final_answer = "1-bromo-4-phenylbutane"
print(f"<<<{final_answer}>>>")