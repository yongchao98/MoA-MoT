def find_smallest_n_for_reducible_space():
    """
    This function determines the smallest integer n for which a reducible n-point
    topological space exists and prints the justification.
    """
    n = 2
    
    # Represent the sets using numbers for clarity
    point1 = 1
    point2 = 2

    # Create string representations for printing
    space_X = f"{{{point1}, {point2}}}"
    subset_Z1 = f"{{{point1}}}"
    subset_Z2 = f"{{{point2}}}"

    print("A topological space is reducible (not irreducible) if it is a union of its proper closed subsets.")
    print("We are looking for the smallest non-negative integer 'n' for which an n-point space can be reducible.\n")
    print("- For n=0 and n=1, any topological space is irreducible.")
    print(f"- For n={n}, we can construct a reducible space.\n")
    
    print(f"Consider the {n}-point space X = {space_X} with the discrete topology.")
    print("In the discrete topology, every subset is a closed set.")
    print(f"Thus, Z1 = {subset_Z1} and Z2 = {subset_Z2} are proper closed subsets of X.")
    
    print("\nThe space X can be expressed as the union of these two proper closed subsets:")
    
    # Print the final equation with each number explicitly shown
    print(f"{space_X} = {subset_Z1} U {subset_Z2}")
    
    print(f"\nSince a 2-point space can be reducible, the smallest such integer is 2.")

find_smallest_n_for_reducible_space()