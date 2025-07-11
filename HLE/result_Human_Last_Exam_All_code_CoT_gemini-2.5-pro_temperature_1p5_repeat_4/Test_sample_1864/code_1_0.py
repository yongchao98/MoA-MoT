import sys

def calculate_flow_controls():
    """
    Calculates the number of essential controls for a flow cytometry sorting experiment.
    
    This calculation includes:
    1. An unstained control for baseline fluorescence.
    2. Single-stain controls for compensation (one for each color).
    3. Fluorescence Minus One (FMO) controls for accurate gating (one for each color).
    """
    
    # The experiment uses five distinct fluorescent channels.
    num_channels = 5
    
    # 1. Unstained Control: Always one sample.
    unstained_controls = 1
    
    # 2. Single-Stain Compensation Controls: One for each channel.
    compensation_controls = num_channels
    
    # 3. Fluorescence Minus One (FMO) Controls: One for each channel, essential for sorting.
    fmo_controls = num_channels
    
    # Calculate the total number of essential control tubes.
    total_controls = unstained_controls + compensation_controls + fmo_controls
    
    print("For a robust 5-channel flow cytometry sorting experiment, you need the following controls:")
    print(f"\n- Unstained Controls: {unstained_controls}")
    print(f"- Single-Stain Compensation Controls: {compensation_controls}")
    print(f"- Fluorescence Minus One (FMO) Controls: {fmo_controls}")
    
    print("\nThis leads to a total number of controls calculated as follows:")
    # The final equation prints each number used in the calculation.
    print(f"{unstained_controls} + {compensation_controls} + {fmo_controls} = {total_controls}")
    
    # In a real script, you might return the value. For this task, we print.
    # We also add the final answer in the required format for the platform.
    # Note: The '<<<...>>>' format is for the platform, not for a typical Python script.
    # We'll print it to standard error to distinguish it from the script's primary output.
    final_answer_for_platform = f"<<<{total_controls}>>>"
    print(final_answer_for_platform, file=sys.stderr)

if __name__ == '__main__':
    calculate_flow_controls()