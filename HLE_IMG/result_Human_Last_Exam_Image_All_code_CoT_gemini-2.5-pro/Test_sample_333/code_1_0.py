def solve_airfoil_matching():
    """
    This function prints the step-by-step reasoning and the final solution
    for the airfoil matching problem.
    """
    
    # The final pairings determined by the analysis
    pairings = {
        'A': 6,
        'B': 2,
        'C': 1,
        'D': 3,
        'E': 8,
        'F': 4,
        'G': 5,
        'H': 7
    }

    print("Analysis complete. The determined pairings are as follows:")
    
    # Print each individual pairing
    for airfoil, plot in pairings.items():
        print(f"Airfoil {airfoil} -> Plot {plot}")

    # Construct the final answer string by concatenating the plot numbers
    final_sequence = "".join(str(pairings[airfoil]) for airfoil in sorted(pairings.keys()))
    
    print("\nThe final sequence of plot numbers for airfoils A-H is:")
    print(final_sequence)

solve_airfoil_matching()