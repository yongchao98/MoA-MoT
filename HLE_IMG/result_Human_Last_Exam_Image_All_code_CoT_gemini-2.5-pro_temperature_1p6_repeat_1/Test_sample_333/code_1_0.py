def solve_airfoil_matching():
    """
    This function explains and prints the solution for matching airfoil shapes to pressure plots.
    """
    print("Based on aerodynamic principles, we can match each airfoil shape (A-H) to its pressure plot (1-8).")
    print("The matching is determined by analyzing camber (related to lift) and thickness/leading edge shape (related to pressure peak).\n")

    # Define the pairings as a dictionary
    pairings = {
        'A': 1,
        'B': 8,
        'C': 4,
        'D': 2,
        'E': 5,
        'F': 6,
        'G': 3,
        'H': 7
    }

    print("The individual pairings are:")
    for airfoil, plot in pairings.items():
        print(f"Airfoil {airfoil} -> Plot {plot}")

    # Generate the final answer string in alphabetical order of airfoils
    final_sequence = "".join(str(pairings[key]) for key in sorted(pairings.keys()))

    print("\nThe final sequence of plot numbers corresponding to airfoils A-H is:")
    print(final_sequence)


if __name__ == "__main__":
    solve_airfoil_matching()