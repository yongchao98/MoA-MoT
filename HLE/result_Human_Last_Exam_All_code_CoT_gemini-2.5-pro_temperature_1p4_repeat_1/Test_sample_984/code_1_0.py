def calculate_efficiency():
    """
    Calculates the efficiency of obtaining paired heavy and light chains
    from the sorted wells in the described experiment.
    """
    # Total number of wells sorted with doublets
    sorted_wells = 95
    
    # Number of wells that yielded paired heavy and light chains
    successful_wells = 4
    
    # The problem mentions a 96-well plate with one well as a negative control.
    # So, the total number of sorted wells is 96 - 1 = 95.
    
    # Calculate the efficiency as a percentage
    efficiency = (successful_wells / sorted_wells) * 100
    
    print(f"Experimental Details:")
    print(f"Total wells sorted: {sorted_wells}")
    print(f"Successful wells with paired chains: {successful_wells}")
    print(f"Calculated Efficiency: {efficiency:.2f}%")
    
    print("\nConclusion:")
    print("The low efficiency of approximately 4.21% is most likely because true antigen-specific B cell:tumor cell interactions are rare events.")
    print("Many of the sorted doublets were likely random, non-specific cell associations rather than true binding pairs, leading to failure in the subsequent PCR step for specific chains.")

calculate_efficiency()