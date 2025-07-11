import sys
# Redirect print to a string to control output format for the final answer
from io import StringIO
original_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

def analyze_protein_folding():
    """
    Analyzes DLS data for MAB13 folding to determine the best conditions
    and evaluates the given answer choices.
    """
    print("Step 1: Define the experimental data and identify the monomer.")
    print("Based on the HEK293 cell data (95% at 7.1 nm), the correctly folded monomer has a radius of 7.1 nm.")
    print("A higher percentage of the 7.1 nm species indicates better protein quality and less aggregation.\n")

    # Store monomer percentages for each condition
    conditions = {
        'E_coli_37C': 0,
        'E_coli_18C': 20,
        'E_coli_18C_HP70': 85,  # Using the better of the two results
        'HEK293_37C': 95,
        'E_coli_37C_GFP': 0,
        'E_coli_18C_MBP': 60
    }

    print("Step 2: Evaluate each answer choice based on the data.\n")

    # --- Evaluation of Choice A ---
    print("--- Evaluating Choice A ---")
    print("Statement: Fusion of another protein to the N-terminal part of MAB13 does not help in the folding process of MAB13.")
    # Check if MBP fusion helped compared to just lowering the temperature
    mbp_improves = conditions['E_coli_18C_MBP'] > conditions['E_coli_18C']
    print(f"Does MBP fusion help? Checking if {conditions['E_coli_18C_MBP']}% > {conditions['E_coli_18C']}%: This is {mbp_improves}.")
    # Check if GFP fusion helped compared to the baseline at the same temperature
    gfp_improves = conditions['E_coli_37C_GFP'] > conditions['E_coli_37C']
    print(f"Does GFP fusion help? Checking if {conditions['E_coli_37C_GFP']}% > {conditions['E_coli_37C']}%: This is {gfp_improves}.")
    # The statement claims fusion does NOT help, but MBP does.
    print("Result: Since MBP fusion improved folding, statement A is FALSE.\n")

    # --- Evaluation of Choice B ---
    print("--- Evaluating Choice B ---")
    print("Statement: Both lower expression temperature and fusion to GFP improve the quality of MAB13.")
    lower_temp_improves = conditions['E_coli_18C'] > conditions['E_coli_37C']
    print(f"Does lower temperature help? Checking if {conditions['E_coli_18C']}% > {conditions['E_coli_37C']}%: This is {lower_temp_improves}.")
    # gfp_improves was already calculated
    print(f"Does GFP fusion help? We already know this is {gfp_improves}.")
    print("Result: Since GFP fusion did not improve quality, statement B is FALSE.\n")

    # --- Evaluation of Choice C ---
    print("--- Evaluating Choice C ---")
    print("Statement: Fusion to MBP improves the folding process of MAB13; MAB13 folds properly in Escherichia coli at 37°C.")
    # mbp_improves was already calculated
    print(f"Does MBP fusion improve folding? We already know this is {mbp_improves}.")
    # Check if it folds properly at 37C in E. coli
    folds_properly_37C = conditions['E_coli_37C'] > 85 # Using a high threshold for "properly"
    print(f"Does MAB13 fold properly in E. coli at 37°C? Checking if {conditions['E_coli_37C']}% is a high value: This is {folds_properly_37C}.")
    print("Result: MAB13 shows 0% proper folding at 37°C in E. coli. Statement C is FALSE.\n")

    # --- Evaluation of Choice D ---
    print("--- Evaluating Choice D ---")
    print("Statement: Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.")
    # First part is the same as A's check
    print("Does adding a fusion improve folding? Yes, MBP fusion did.")
    # Second part is the same as C's check
    print("Can MAB13 fold properly at 37°C? No, not in E. coli (0% monomer).")
    print("Result: Statement D is FALSE.\n")

    # --- Evaluation of Choice E ---
    print("--- Evaluating Choice E ---")
    print("Statement: Both GFP and HP70 do not facilitate the folding of MAB13.")
    # gfp_improves is False, so "GFP does not facilitate" is True
    print(f"Does GFP facilitate folding? No, the monomer percent did not increase from {conditions['E_coli_37C']}%.")
    hp70_improves = conditions['E_coli_18C_HP70'] > conditions['E_coli_18C']
    print(f"Does HP70 facilitate folding? Checking if {conditions['E_coli_18C_HP70']}% > {conditions['E_coli_18C']}%: This is {hp70_improves}.")
    print("Result: Since HP70 does facilitate folding, statement E is FALSE.\n")

    # --- Evaluation of Choice F ---
    print("--- Evaluating Choice F ---")
    print("Statement: HP70 facilitates the folding process of MAB13 at 18°C and 37°C, MBP and lower temperature improve the folding process of MAB13.")
    # Check each part of the statement
    # hp70_improves was already calculated
    print(f"Does HP70 facilitate folding at 18°C? Yes, monomer percent increased from {conditions['E_coli_18C']}% to {conditions['E_coli_18C_HP70']}%.")
    print("Does HP70 facilitate folding at 37°C? This data is not provided.")
    # mbp_improves was already calculated
    print(f"Does MBP improve the folding process? Yes, monomer percent increased from {conditions['E_coli_18C']}% to {conditions['E_coli_18C_MBP']}%.")
    # lower_temp_improves was already calculated
    print(f"Does lower temperature improve the folding process? Yes, monomer percent increased from {conditions['E_coli_37C']}% to {conditions['E_coli_18C']}%.")
    print("Result: All testable parts of this statement are TRUE. While the data for HP70 at 37°C is missing, all other statements (A-E) are demonstrably false. Therefore, F is the best possible answer.\n")

    print("Step 3: Final Conclusion.")
    print("Based on the analysis, statement F is the most accurate description of the results.")
    
analyze_protein_folding()

# Restore original stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())
print("<<<F>>>")