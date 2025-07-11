import sys
import io

# Helper function to redirect stdout to capture the output for final formatting
def get_script_output(script_func):
    old_stdout = sys.stdout
    sys.stdout = new_stdout = io.StringIO()
    script_func()
    sys.stdout = old_stdout
    return new_stdout.getvalue()

def solve():
    # --- Experiment 1 Data ---
    # Plaque-forming units (cfu) per microliter
    # Bacteria without RP system
    cfu_wt_noRP = 100000
    cfu_deltaXY_noRP = 100000
    # Bacteria with RP system
    cfu_wt_RP = 80000
    cfu_deltaXY_RP = 40000

    # --- Experiment 2 Data ---
    # Detection of 500 Da molecule 60 minutes post-infection
    # True = detected, False = not detected
    # Sample 1: Bacteria with RP + PhageDE3-wt (has XY)
    s1_molecule_detected = True
    # Sample 2: Bacteria with RP + PhageDE3-deltaXY (no XY)
    s2_molecule_detected = False
    # Sample 3: Bacteria without RP + PhageDE3-wt (has XY)
    s3_molecule_detected = False
    # Sample 4: Bacteria without RP + PhageDE3-deltaXY (no XY)
    s4_molecule_detected = False

    # --- Analysis ---

    print("Step 1: Analyzing the effect of the RP system on bacterial resistance.")
    print(f"To assess resistance, we compare the phage's success against bacteria with and without the RP system.")
    print(f"For the phage without operon XY, the CFU dropped from {cfu_deltaXY_noRP} (without RP) to {cfu_deltaXY_RP} (with RP).")
    print("Conclusion 1: Since the phage is less effective against bacteria with the RP system, the RP system increases bacterial resistance.\n")

    print("Step 2: Analyzing the phage's counter-defense (Operon XY).")
    print("In the presence of the RP system, we compare the phage with and without the operon XY.")
    print(f"Phage-wt (with XY) CFU = {cfu_wt_RP}")
    print(f"Phage-deltaXY (without XY) CFU = {cfu_deltaXY_RP}")
    print(f"Conclusion 2: The phage with operon XY ({cfu_wt_RP} cfu) is more effective than the phage without it ({cfu_deltaXY_RP} cfu) when the RP system is present. This suggests operon XY is a counter-defense system.\n")
    
    print("Step 3: Analyzing the production of the 500 Da molecule.")
    print("The molecule with mass 500 Da was detected ONLY under the following conditions:")
    print(f"  - Bacteria with RP system: YES")
    print(f"  - Phage with operon XY: YES")
    print("The molecule was NOT detected if either the RP system or the XY operon was absent.")
    print("Conclusion 3: The enzymes from operon XY synthesize the 500 Da molecule, and this process requires a component from the bacterial RP system.\n")
    
    print("Step 4: Final Evaluation")
    print("Based on the analysis, the full story is:")
    print("  1. The RP system provides resistance to the phage (Conclusion 1).")
    print("  2. The phage's XY operon counters this resistance (Conclusion 2).")
    print("  3. The mechanism of the counter-defense involves the XY enzymes creating a 500 Da molecule, a process which depends on the presence of the RP system (Conclusion 3).")
    
    print("\nStatement H ('System RP increases the resistance of the bacteria against phageDE3 because the enzymes XY1 or XY2 can synthesize their products only in the presence of the system RP.') best summarizes these interconnected findings. It correctly identifies that RP increases resistance and correctly states the condition required for the synthesis of the phage's counter-defense product.")

# Run the analysis and capture output
output = get_script_output(solve)
print(output)
print("<<<H>>>")