def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a multi-color flow cytometry
    sorting experiment based on user-provided parameters.
    """
    # --- Experiment Parameters ---
    channels = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_channels = len(channels)
    # The mention of "Streptavidin beads" implies a biotin-streptavidin system is
    # used for detection, which requires a specific control.
    uses_streptavidin = True

    # --- Control Calculation ---
    # 1. Unstained control: a baseline for autofluorescence.
    unstained_control_count = 1

    # 2. Single-stain controls: for compensation of spectral overlap. One for each color.
    single_stain_control_count = num_channels

    # 3. Fluorescence Minus One (FMO) controls: for accurate gating, essential for sorting. One for each color.
    fmo_control_count = num_channels

    # 4. Streptavidin reagent control: checks for non-specific binding of the streptavidin conjugate.
    streptavidin_control_count = 1 if uses_streptavidin else 0

    # --- Total Calculation ---
    total_controls = (unstained_control_count +
                      single_stain_control_count +
                      fmo_control_count +
                      streptavidin_control_count)

    # --- Output Explanation ---
    print("For a robust five-channel flow cytometry sorting experiment, you should prepare the following essential controls:")
    print("-" * 80)

    print(f"1. Unstained Control:")
    print(f"   - Purpose: To measure the baseline autofluorescence of your sample.")
    print(f"   - Number of tubes: {unstained_control_count}\n")

    print(f"2. Single-Stain Compensation Controls:")
    print(f"   - Purpose: To calculate the compensation matrix and correct for spectral spillover.")
    print(f"   - One control is needed for each fluorochrome: {', '.join(channels)}.")
    print(f"   - Number of tubes: {single_stain_control_count}\n")

    print(f"3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Purpose: To accurately set positive gates by visualizing signal spread from all other")
    print(f"     fluorochromes. This is critical for sorting experiments.")
    print(f"   - One FMO is prepared for each color in your panel.")
    print(f"   - Number of tubes: {fmo_control_count}\n")

    if uses_streptavidin:
        print(f"4. Streptavidin Reagent Control:")
        print(f"   - Purpose: To control for non-specific binding of the fluorescent streptavidin conjugate.")
        print(f"     This is essential because a biotin-streptavidin system requires this check.")
        print(f"   - Number of tubes: {streptavidin_control_count}\n")

    print("-" * 80)
    print("Total Number of Essential Controls:")
    
    # Final equation with each number explicitly shown
    equation = (f"{unstained_control_count} (Unstained) + "
                f"{single_stain_control_count} (Single-Stain) + "
                f"{fmo_control_count} (FMO) + "
                f"{streptavidin_control_count} (Streptavidin) = {total_controls}")
    print(equation)
    print("-" * 80)

if __name__ == '__main__':
    calculate_flow_controls()