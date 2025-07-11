def get_ideal_ni_ce_ratio():
    """
    Summarizes the ideal Ni/Ce ratio for Ni-Ceria catalysts based on scientific literature.
    """
    
    # Define the optimal range and example values based on literature consensus
    lower_bound_ratio = 0.1
    upper_bound_ratio = 0.3
    example_ni_part = 1
    example_ce_part = 9 # Corresponds to a Ni/Ce ratio of 1/9 or ~0.11

    print("--- Ideal Ni/Ce Ratio for Ni-Ceria Catalysts ---")
    print("\nThe ideal Ni/Ce ratio is a range designed to maximize the active sites at the Nickel-Ceria interface.")
    
    print("\n--- Optimal Range ---")
    print(f"Based on extensive research for both Water-Gas Shift (WGS) and Water Splitting (WS) reactions, the optimal Ni/Ce atomic ratio generally falls within the range of:")
    print(f"Lower Bound: {lower_bound_ratio}")
    print(f"Upper Bound: {upper_bound_ratio}")
    print(f"\nThis can be expressed in various ways. For example, a common high-performance composition is a molar ratio of approximately {example_ni_part} part Ni to {example_ce_part} parts Ce.")

    print("\n--- Scientific Rationale ---")
    print("1. Maximizing the Interface: The synergy between Nickel and Ceria is crucial. The reaction primarily occurs at the interface between the Ni nanoparticles and the Ceria (CeO2) support. A ratio within the optimal range ensures a high density of these active sites.")
    print("2. High Ni Dispersion: Lower Ni content (e.g., Ni/Ce < 0.3) promotes the formation of small, highly dispersed Ni nanoparticles on the Ceria support. This maximizes the surface area of the active metal.")
    print("3. Preventing Sintering: At higher Ni concentrations (e.g., Ni/Ce > 0.3), the Ni nanoparticles tend to agglomerate or 'sinter' into larger particles during the reaction. This reduces the active surface area and the number of beneficial interface sites, leading to a drop in catalytic performance.")
    
    print("\n--- Conclusion ---")
    print("The key is to achieve a balance: enough Ni to provide catalytic sites, but not so much that it agglomerates. The synthesis method plays a critical role in achieving this ideal nanostructure.")
    
# Execute the function to print the information
get_ideal_ni_ce_ratio()