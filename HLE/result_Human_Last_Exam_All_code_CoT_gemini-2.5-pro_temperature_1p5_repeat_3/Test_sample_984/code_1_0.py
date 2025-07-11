def calculate_efficiency():
    """
    Calculates and prints the efficiency of the cell sorting experiment.
    """
    # Total number of wells into which doublets were sorted
    total_sorted_wells = 95
    
    # Number of wells that yielded paired heavy and light chains
    successful_wells = 4
    
    # Calculate the efficiency as a percentage
    efficiency = (successful_wells / total_sorted_wells) * 100
    
    print("The experimental efficiency is calculated as follows:")
    # The final equation is printed showing each number
    print(f"({successful_wells} successful wells / {total_sorted_wells} sorted wells) * 100 = {efficiency:.2f}%")

calculate_efficiency()