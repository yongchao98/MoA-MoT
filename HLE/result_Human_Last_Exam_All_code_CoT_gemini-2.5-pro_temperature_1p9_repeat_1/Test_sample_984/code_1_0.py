import sys

def calculate_efficiency():
    """
    Calculates and prints the efficiency of the single-cell sorting experiment.
    """
    # Total number of wells sorted with potential doublets.
    # One of the 96 wells was a negative control.
    total_sorted_wells = 95
    
    # Number of wells that successfully yielded paired heavy and light chains.
    successful_outcomes = 4
    
    # Check for division by zero, although not possible with these hardcoded numbers.
    if total_sorted_wells == 0:
        print("Error: Cannot calculate efficiency with zero sorted wells.", file=sys.stderr)
        return
        
    # Calculate the efficiency.
    efficiency_percentage = (successful_outcomes / total_sorted_wells) * 100
    
    print("Analysis of Experimental Efficiency:")
    print(f"Total number of sorted doublets: {total_sorted_wells}")
    print(f"Number of successful paired chain amplifications: {successful_outcomes}")
    print(f"The efficiency calculation is: ({successful_outcomes} / {total_sorted_wells}) * 100")
    print(f"The resulting efficiency is only {efficiency_percentage:.2f}%.")
    print("\nThis very low efficiency rate suggests that the vast majority of sorted 'doublets' did not represent the specific biological interaction of interest.")

calculate_efficiency()