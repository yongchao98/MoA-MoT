import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

def analyze_probe_solubility():
    """
    Analyzes and compares the estimated solubility of a chemical probe before and after modification.

    This function uses a simplified model based on hydrophilicity scores for different
    molecular fragments to demonstrate why PEGylation is a good strategy to solve
    precipitation issues. A lower total score indicates higher predicted water solubility.
    """

    # These are illustrative scores for a simplified model.
    # Positive scores = hydrophobic (bad for water solubility)
    # Negative scores = hydrophilic (good for water solubility)
    hydrophilicity_scores = {
        "Thioxanthenone_Core": 15.0,  # A large, planar, aromatic system is very hydrophobic.
        "Chlorohexyl_Tail": 8.0,      # A 6-carbon alkyl chain is hydrophobic.
        "Existing_Ethers": -3.0,      # The two existing ether groups provide some hydrophilicity.
        "Amide_Group": -2.0,          # The amide group is hydrophilic, but less so than a PEG chain.
        "PEG_Linker_Swap": -7.0       # Replacing the amide with a short PEG linker provides significant hydrophilicity.
    }

    # --- Part 1: Calculate the score for the Original Probe ---
    core = hydrophilicity_scores["Thioxanthenone_Core"]
    tail = hydrophilicity_scores["Chlorohexyl_Tail"]
    ethers = hydrophilicity_scores["Existing_Ethers"]
    amide = hydrophilicity_scores["Amide_Group"]

    original_probe_score = core + tail + ethers + amide

    print("--- Analysis of Probe Solubility ---")
    print("A lower score corresponds to higher predicted water solubility.\n")

    print("Original Probe Calculation:")
    print("The final score is the sum of the scores of its parts:")
    print(f"Equation: {core} (Core) + {tail} (Tail) + {ethers} (Ethers) + {amide} (Amide)")
    print(f"Final Score for Original Probe = {original_probe_score}\n")

    # --- Part 2: Calculate the score for the Modified Probe ---
    # The modified probe swaps the amide group for a more powerful hydrophilic PEG linker.
    peg_linker = hydrophilicity_scores["PEG_Linker_Swap"]

    modified_probe_score = core + tail + ethers + peg_linker

    print("Modified Probe Calculation (Amide replaced with PEG Linker):")
    print("The final score is the sum of the scores of its new parts:")
    print(f"Equation: {core} (Core) + {tail} (Tail) + {ethers} (Ethers) + {peg_linker} (PEG Linker)")
    print(f"Final Score for Modified Probe = {modified_probe_score}\n")

    # --- Part 3: Conclusion ---
    print("--- Conclusion ---")
    if modified_probe_score < original_probe_score:
        print("Yes, changing the amide group to a PEG group is an excellent strategy.")
        print("The model shows a decrease in the hydrophobicity score from "
              f"{original_probe_score} to {modified_probe_score}.")
        print("This modification increases the overall hydrophilicity of the probe, "
              "which should significantly improve its solubility in the aqueous cell medium and solve the precipitation problem.")
    else:
        # This case is unlikely given the model's setup
        print("The model suggests the modification may not improve solubility. "
              "However, in practice, PEGylation is a proven method for increasing solubility.")

# Run the analysis
analyze_probe_solubility()

# --- Final Output ---
# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = string_buffer.getvalue()
# Print the output
print(output)
# Also return the final conclusion as the answer
print("<<<Yes>>>")