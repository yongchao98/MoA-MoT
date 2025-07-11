def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry experiment.

    This script determines the minimum number of controls required for instrument setup
    and spectral compensation for a 5-color sorting experiment using beads.
    """

    # Number of different fluorochromes being used in the experiment.
    num_fluorochromes = 5
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unlabeled Beads Control:
    # This control contains only the beads and is used to set the baseline fluorescence (autofluorescence).
    unlabeled_beads_control = 1

    # 2. Single-Stain Compensation Controls:
    # One control is needed for each fluorochrome to calculate the compensation matrix
    # and correct for spectral overlap.
    single_stain_controls = num_fluorochromes

    # Calculate the total number of essential controls.
    total_essential_controls = unlabeled_beads_control + single_stain_controls

    # Print the explanation and the final result.
    print("Essential Controls for a 5-Color Streptavidin Bead Sorting Experiment:")
    print("-" * 65)
    print(f"1. Unlabeled Beads Control: A tube with just the beads.")
    print(f"   - Purpose: To measure bead autofluorescence and set detector voltages.")
    print(f"   - Count: {unlabeled_beads_control}\n")

    print(f"2. Single-Stain Compensation Controls: One for each fluorochrome.")
    print(f"   - Fluorochromes: {', '.join(fluorochromes)}")
    print(f"   - Purpose: To calculate the compensation matrix for spectral spillover.")
    print(f"   - Count: {single_stain_controls}\n")

    print("Total Number of Essential Controls Calculation:")
    print(f"This is the sum of the unlabeled control and the single-stain controls.")
    print(f"Equation: {unlabeled_beads_control} (unlabeled) + {single_stain_controls} (single-stains) = {total_essential_controls}")
    print("-" * 65)
    print(f"You should prepare a minimum of {total_essential_controls} essential controls.")
    print("\nNote: For accurate gating, especially with dim or overlapping populations,")
    print("it is also highly recommended to prepare Fluorescence Minus One (FMO) controls.")


# Run the calculation.
calculate_flow_controls()
<<<6>>>