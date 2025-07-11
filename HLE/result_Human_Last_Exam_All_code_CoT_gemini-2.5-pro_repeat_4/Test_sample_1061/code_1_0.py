import sys

def model_tree_ring_isotope(base_ratio, environmental_factor, sensitivity):
    """
    A simplified model to demonstrate how an environmental factor affects the 13C ratio.
    In this case, the environmental factor is monsoon strength (water availability).
    A stronger monsoon leads to a lower (more negative) 13C ratio.
    
    The equation is: final_ratio = base_ratio - (environmental_factor * sensitivity)
    
    Args:
        base_ratio (float): The starting delta 13C value (in per mil, ‰).
        environmental_factor (float): A normalized value for the strength of the factor.
        sensitivity (float): How much the ratio changes per unit of the factor.
        
    Returns:
        float: The resulting 13C ratio.
    """
    return base_ratio - (environmental_factor * sensitivity)

def main():
    """
    Main function to run the simulation and explain the answer.
    """
    print("Analyzing the predominant factor for declining 13C ratio in tree rings (1886-1990 AD).")
    print("The key principle is that greater water availability (e.g., from a stronger monsoon) allows trees to discriminate more against the heavy 13C isotope, resulting in a lower 13C ratio in their rings.")
    print("-" * 80)

    # --- Model Parameters ---
    # A hypothetical base 13C ratio (in per mil, ‰) under neutral conditions.
    base_c13_ratio = -25.0
    # A sensitivity factor: how much the 13C ratio changes per unit of monsoon strength.
    sensitivity_factor = 0.04

    # --- Simulation for the start of the period (1886) ---
    year_start = 1886
    # Let's assume a baseline monsoon strength at the start.
    monsoon_strength_start = 10.0
    
    # Calculate the 13C ratio for the start year
    ratio_start = model_tree_ring_isotope(base_c13_ratio, monsoon_strength_start, sensitivity_factor)
    
    print(f"Simulation for Year: {year_start}")
    print("Equation: final_ratio = base_ratio - (monsoon_strength * sensitivity)")
    # Using sys.stdout.write to prevent the print function from adding extra spaces
    sys.stdout.write(f"Result: {ratio_start:.3f} = {base_c13_ratio} - ({monsoon_strength_start} * {sensitivity_factor})\n")
    print("-" * 80)
    
    # --- Simulation for the end of the period (1990) ---
    year_end = 1990
    # Assume a trend of a strengthening monsoon over the century.
    monsoon_strength_end = 45.0
    
    # Calculate the 13C ratio for the end year
    ratio_end = model_tree_ring_isotope(base_c13_ratio, monsoon_strength_end, sensitivity_factor)
    
    print(f"Simulation for Year: {year_end}")
    print("Equation: final_ratio = base_ratio - (monsoon_strength * sensitivity)")
    sys.stdout.write(f"Result: {ratio_end:.3f} = {base_c13_ratio} - ({monsoon_strength_end} * {sensitivity_factor})\n")
    print("-" * 80)
    
    print(f"Conclusion: As the modeled monsoon strength increased from {monsoon_strength_start} to {monsoon_strength_end},")
    print(f"the 13C ratio declined from {ratio_start:.3f}‰ to {ratio_end:.3f}‰.")
    print("This demonstrates that changes in the SE Asia monsoon are a plausible predominant factor for the observed trend.")

if __name__ == "__main__":
    main()

<<<E>>>