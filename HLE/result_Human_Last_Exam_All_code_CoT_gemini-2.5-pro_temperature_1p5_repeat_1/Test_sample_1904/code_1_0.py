import itertools

def metric(p1, p2):
    """
    Defines the discrete ultrametric on a set of points.
    d(x, y) = 1 if x != y, and 0 if x == y.
    """
    if p1 == p2:
        return 0
    else:
        return 1

def calculate_diameter(subset):
    """
    Calculates the diameter of a given subset of points.
    The diameter is the supremum (or maximum in this finite case)
    of distances between all pairs of points in the subset.
    """
    if not subset:
        return 0 # Diameter of an empty set is 0
    if len(subset) == 1:
        return 0

    max_dist = 0
    # Check all pairs of points in the subset
    for p1, p2 in itertools.combinations(subset, 2):
        dist = metric(p1, p2)
        if dist > max_dist:
            max_dist = dist
    return max_dist

def main():
    """
    Main function to demonstrate the calculation of diameters for our constructed space.
    """
    # Let's model the infinite space X with a finite sample of 10 points.
    X = set(range(10))
    print(f"Sample space X = {X}\n")

    # Define some non-empty closed subsets. In a discrete space, any subset is closed.
    subsets_to_test = [
        {5},                  # A singleton set
        {3, 8},               # A set with two points
        {0, 1, 2},            # A set with three points
        {0, 2, 4, 6, 8},      # A larger set
        X                     # The entire sample space
    ]

    print("Calculating diameters of various non-empty closed subsets...")
    
    unique_diameters = set()

    for i, subset in enumerate(subsets_to_test):
        diam = calculate_diameter(subset)
        unique_diameters.add(diam)
        print(f"  - Subset A_{i+1} = {subset}:")
        if len(subset) > 1:
            # Show an example calculation
            p1, p2 = next(itertools.combinations(subset, 2))
            print(f"    Example distance: d({p1}, {p2}) = {metric(p1, p2)}")
        print(f"    Diameter = {diam}\n")

    print("Summary of results:")
    print(f"The set of unique diameters found is: {unique_diameters}")
    
    num_components = len(unique_diameters)
    print("The number of possible diameters determines the number of connected components.")
    # The final equation part requested by the prompt
    print(f"Smallest number of connected components = |Set of Diameters| = {num_components}")


if __name__ == "__main__":
    main()