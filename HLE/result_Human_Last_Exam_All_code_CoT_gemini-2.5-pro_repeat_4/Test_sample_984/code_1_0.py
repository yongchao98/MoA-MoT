def calculate_efficiency():
    """
    This function calculates the efficiency of the cell sorting experiment
    based on the numbers provided in the problem description.
    """
    # Number of wells where doublets were sorted into
    total_sorted_wells = 95
    
    # Number of wells that yielded paired heavy and light chains
    successful_wells = 4
    
    # Calculate the percentage efficiency
    efficiency = (successful_wells / total_sorted_wells) * 100
    
    print("The experimental efficiency is calculated as follows:")
    # The prompt requires printing each number in the final equation.
    print(f"Equation: ({successful_wells} successful wells / {total_sorted_wells} total sorted wells) * 100")
    print(f"Resulting Efficiency: {efficiency:.2f}%")
    
    print("\nThis low efficiency highlights that the event of interest was rare.")
    print("The most likely explanation is that true, antigen-specific interactions are rare, and most of the sorted doublets were random cell pairings rather than specific binding events.")

calculate_efficiency()