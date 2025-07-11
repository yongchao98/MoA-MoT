import math

def calculate_and_print_shannon(site_name, counts_dict):
    """
    Calculates and prints the Shannon diversity calculation for a given site.
    
    Args:
        site_name (str): The name of the site.
        counts_dict (dict): A dictionary of species names and their counts.
        
    Returns:
        float: The calculated Shannon diversity index.
    """
    counts = list(counts_dict.values())
    N = sum(counts)
    
    print(f"--- Calculation for {site_name} ---")
    print(f"True Species Counts: {counts_dict}")
    print(f"Total individuals (N) = {N}")

    # If there is 0 or 1 species, diversity is 0
    if len(counts) <= 1:
        shannon_index = 0.0
        # The formula for one species simplifies to 0
        p = 1.0
        print(f"Shannon's H = - ( {p:.2f} * ln({p:.2f}) ) = {shannon_index:.2f}")
        print("-" * 35 + "\n")
        return shannon_index

    proportions = [count / N for count in counts]
    
    # Build the equation string with all numbers
    equation_parts = []
    for p in proportions:
        # We use .4f for precision in the equation parts
        equation_parts.append(f"{p:.4f} * ln({p:.4f})")
    equation_str = f"- ( {' + '.join(equation_parts)} )"
    
    # Calculate the final value
    shannon_index = -sum(p * math.log(p) for p in proportions if p > 0)
    
    print(f"Shannon's H = {equation_str}")
    print(f"Shannon's H = {shannon_index:.2f}")
    print("-" * 35 + "\n")
    
    return shannon_index

def main():
    """
    Main function to consolidate data and run calculations.
    """
    # Consolidated abundance data based on true species identification
    site_data = {
        'Site 1': {'Harmonia axyridis': 4},
        'Site 2': {'Harmonia axyridis': 23},
        'Site 3': {'Harmonia axyridis': 4, 'Species F': 4}
    }
    
    # Calculate diversity for each site
    diversities = []
    for site, counts in site_data.items():
        diversity = calculate_and_print_shannon(site, counts)
        diversities.append(diversity)
        
    # Format and print the final answer as requested
    final_answer = ", ".join([f"{d:.2f}" for d in diversities])
    print("The insect diversity of each site to two decimal places, separated by \", \" is:")
    print(final_answer)

if __name__ == "__main__":
    main()