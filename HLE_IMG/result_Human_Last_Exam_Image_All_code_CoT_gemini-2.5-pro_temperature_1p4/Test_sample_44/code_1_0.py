import sys
import io

# Redirect stdout to capture the print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_nmr_puzzle():
    """
    Analyzes 1H NMR data and molecular structures to identify the correct compound.
    """
    # Step 1: Analyze the provided 1H NMR data.
    print("### Step 1: Analyzing the 1H NMR Data ###\n")
    nmr_integrations = [1, 1, 1, 1, 1, 1, 1, 4, 4, 2, 2, 2]
    nmr_shifts = [8.19, 7.79, 7.47, 7.38, 6.98, 6.63, 6.61, 4.19, 3.63, 3.21, 2.83, 1.98]
    
    nmr_total_h = sum(nmr_integrations)
    nmr_aromatic_h = sum(h for i, h in enumerate(nmr_integrations) if nmr_shifts[i] > 6.0)
    nmr_aliphatic_h = sum(h for i, h in enumerate(nmr_integrations) if nmr_shifts[i] < 5.0)

    print("From the provided NMR spectrum:")
    print(f" - Total integration = {nmr_total_h} protons")
    print(f" - Aromatic protons (in the region > 6.0 ppm) = {nmr_aromatic_h}H")
    print(f" - Aliphatic protons (in the region < 5.0 ppm) = {nmr_aliphatic_h}H")
    print("\n" + "-"*40 + "\n")

    # Step 2: Analyze the proton counts for each candidate structure.
    print("### Step 2: Analyzing the Candidate Structures ###\n")
    # Structure A: Pyridyl-piperazine + Tetrahydroquinoline-thiosemicarbazone
    A_aromatic = 7  # 3 from quinoline, 4 from pyridine
    A_aliphatic = 14 # 6 from quinoline, 8 from piperazine
    A_nh = 1
    A_total = A_aromatic + A_aliphatic + A_nh
    print(f"Compound A: Total H = {A_total}, Aromatic H = {A_aromatic}, Aliphatic+NH H = {A_aliphatic + A_nh}")

    # Structure C: Phenyl-piperazine + Tetrahydroquinoline-thiosemicarbazone
    C_aromatic = 8  # 3 from quinoline, 5 from phenyl
    C_aliphatic = 14 # 6 from quinoline, 8 from piperazine
    C_nh = 1
    C_total = C_aromatic + C_aliphatic + C_nh
    print(f"Compound C: Total H = {C_total}, Aromatic H = {C_aromatic}, Aliphatic+NH H = {C_aliphatic + C_nh}")

    # Structures B, D, E are complexes with two ligands
    BDE_total = A_total * 2 # at least 44 protons
    print(f"Compounds B, D, E: These are metal complexes containing two ligands, resulting in a much higher proton count (>= {BDE_total}H). They are inconsistent with the NMR data's total of {nmr_total_h}H.")
    print("\n" + "-"*40 + "\n")
    
    # Step 3: Compare data and structures, accounting for potential missing NH proton
    print("### Step 3: Comparison and Conclusion ###\n")
    print(f"The NMR data shows {nmr_total_h} total protons, while Compound A has {A_total} and C has {C_total}.")
    print("A common reason for a 1H discrepancy is the exchange of an acidic NH proton with solvent, making it 'invisible' to the NMR.")
    print("Let's re-evaluate assuming the NH proton is not observed:\n")

    A_total_no_nh = A_total - A_nh
    
    print(f"Hypothesis: The data is from Compound A, with the NH proton not observed.")
    print(f" - Expected Total H = {A_total} - 1 (NH) = {A_total_no_nh}")
    print(f" - Expected Aromatic H = {A_aromatic}")
    print(f" - Expected Aliphatic H = {A_aliphatic}")
    
    print("\nComparing this hypothesis with the NMR data:")
    print(f"                      Structure A (no NH)   |   NMR Data")
    print(f"Total H:              {A_total_no_nh:<23} |   {nmr_total_h}")
    print(f"Aromatic H:           {A_aromatic:<23} |   {nmr_aromatic_h}")
    print(f"Aliphatic H:          {A_aliphatic:<23} |   {nmr_aliphatic_h}")

    # Final check for matching
    if A_total_no_nh == nmr_total_h and A_aromatic == nmr_aromatic_h and A_aliphatic == nmr_aliphatic_h:
        print("\nThe counts match perfectly. The NMR data corresponds to structure A.")
        # The answer choices are A. B, B. D, C. E, D. C, E. A.
        # The compound is A, which is option E.
        final_answer = "<<<E>>>"
    else:
        final_answer = "<<<No Match Found>>>"
        
    return final_answer

# Run the solver and get the final answer
final_answer_string = solve_nmr_puzzle()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print(final_answer_string)
