def find_ideal_ratio():
    """
    This function explains and provides the ideal Ni/Ce ratio for catalytic reactions
    based on findings from scientific literature.
    """

    print("The ideal Ni/Ce ratio for Ni-Ceria catalysts is not a single, universal value.")
    print("It depends on factors like synthesis method, particle size, and specific reaction conditions (e.g., temperature, pressure).")
    print("\nHowever, extensive research has identified an optimal range that maximizes catalytic performance for reactions like Water Gas Shift (WGS) and water splitting.")
    print("\nThe goal is to achieve high dispersion of Ni atoms on the Ceria (CeO₂) support, maximizing the active Ni-CeO₂ interface.")
    
    # Define the optimal range based on literature
    # This is typically expressed as an atomic or molar ratio.
    lower_bound_ratio = 0.1
    upper_bound_ratio = 0.2
    
    print("\n--- Key Finding ---")
    print("A low Ni loading is generally preferred to prevent the agglomeration (sintering) of Ni particles and to promote strong metal-support interactions.")
    print(f"The most frequently reported optimal atomic molar ratio of Ni/Ce is in the range of {lower_bound_ratio} to {upper_bound_ratio}.")
    print(f"This corresponds to an atomic composition of 1 Ni atom for every 5 to 10 Ce atoms.")
    
    # Example of a highly cited ratio within this range
    example_ratio = 0.11
    
    print(f"\nFor example, a ratio of approximately {example_ratio} (often prepared as 10 wt% Ni on CeO₂) is frequently cited for achieving high catalytic activity and stability.")

# Execute the function to print the information
find_ideal_ratio()