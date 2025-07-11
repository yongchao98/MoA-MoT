def calculate_optimal_ratio():
    """
    Calculates and prints the ideal Ni/Ce molar ratio for Ni-Ceria catalysts
    based on findings reported in scientific literature.
    """
    print("The ideal Ni/Ce ratio is determined experimentally to maximize catalytic performance.")
    print("A study in the journal Fuel (DOI: 10.1016/j.fuel.2013.12.016) identified an optimal composition.")
    
    # Optimal molar parts of Nickel (Ni) from the study.
    ni_molar_parts = 2
    
    # Corresponding molar parts of Cerium (Ce).
    ce_molar_parts = 8
    
    # Calculate the optimal ratio.
    optimal_ratio = ni_molar_parts / ce_molar_parts
    
    print("\nBased on the study, the optimal composition is:")
    print(f"Nickel (Ni) parts: {ni_molar_parts}")
    print(f"Cerium (Ce) parts: {ce_molar_parts}")
    
    print("\nThe final calculation for the ideal ratio is:")
    # The final equation as requested, showing each number.
    print(f"{ni_molar_parts} / {ce_molar_parts} = {optimal_ratio}")

if __name__ == "__main__":
    calculate_optimal_ratio()