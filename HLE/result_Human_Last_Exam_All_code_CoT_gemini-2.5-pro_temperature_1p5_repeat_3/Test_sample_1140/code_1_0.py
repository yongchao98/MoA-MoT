import math

def calculate_synthesis_yields():
    """
    Calculates and compares the theoretical yields for direct SPPS vs. an NCL fragment.

    This script illustrates why Native Chemical Ligation (NCL) is a superior strategy
    for synthesizing long peptides by showing the dramatic improvement in yield when
    synthesizing smaller fragments.
    """
    # --- Parameters ---
    total_peptide_length = 100  # Total amino acids in the final peptide
    coupling_efficiency = 0.99  # A typical per-step efficiency for SPPS
    num_fragments_ncl = 2     # We plan to use NCL with 2 fragments

    # --- Calculations ---

    # 1. Direct Synthesis of the full 100aa peptide via SPPS
    # The number of coupling steps is the length - 1.
    direct_spps_steps = total_peptide_length - 1
    direct_spps_yield = math.pow(coupling_efficiency, direct_spps_steps)

    # 2. Synthesis of a single fragment for NCL
    # We split the 100aa peptide into 2 smaller, more manageable fragments.
    fragment_length = total_peptide_length // num_fragments_ncl
    fragment_spps_steps = fragment_length - 1
    fragment_spps_yield = math.pow(coupling_efficiency, fragment_spps_steps)

    # --- Output ---
    print("--- Comparing Peptide Synthesis Yields ---\n")
    print(f"Assumptions:")
    print(f"  - Total Peptide Length: {total_peptide_length} aa")
    print(f"  - Per-Step Coupling Efficiency: {coupling_efficiency * 100}%\n")

    print("Method 1: Direct Solid-Phase Peptide Synthesis (SPPS)")
    print("This method synthesizes the full peptide in one continuous run.")
    # The final print statement fulfills the requirement to show each number in the equation.
    print(f"Calculation: {coupling_efficiency} ^ ({total_peptide_length} - 1)")
    print(f"Theoretical Yield: {direct_spps_yield:.3%}\n")

    print(f"Method 2: Native Chemical Ligation (NCL) Approach")
    print(f"This method first synthesizes a smaller {fragment_length}aa fragment.")
    # The final print statement fulfills the requirement to show each number in the equation.
    print(f"Calculation: {coupling_efficiency} ^ ({fragment_length} - 1)")
    print(f"Theoretical Yield for one fragment: {fragment_spps_yield:.3%}\n")
    
    print("--- Conclusion ---")
    print("The yield for a 50aa fragment is significantly higher than for the full 100aa peptide.")
    print("This higher yield makes purification of the fragments and the subsequent ligation step to create the final product much more feasible.")

# Run the calculation
calculate_synthesis_yields()