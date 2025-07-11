import sys

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for the given Calabi-Yau Link.
    """
    # Step 1: Define the weights for the variables (z1, z2, z3, z4, z5)
    weights = [22, 29, 49, 50, 75]

    # Step 2: Determine the quasi-homogeneous degree 'd'.
    # Based on the analysis, most terms have a weighted degree of 225.
    # We assume the polynomial is quasi-homogeneous with this degree.
    # For example, from the term z5^3: d = 3 * 75 = 225.
    d = 225

    # Step 3: Calculate the sum of the weights.
    sum_of_weights = sum(weights)

    # A Calabi-Yau hypersurface is defined by the condition d = sum(weights).
    # Here, 225 == 225, so the condition holds.

    # Step 4: Calculate the Crawley-Nordström invariant.
    # The formula is (d - sum(weights)) / d.
    # Using float division to ensure the result is a float, although it's 0 here.
    if d == 0:
      print("Error: Degree 'd' cannot be zero.", file=sys.stderr)
      return
    invariant = (d - sum_of_weights) / float(d)

    # Print the explanation and the final equation as requested.
    weights_str = " + ".join(map(str, weights))
    print(f"The Crawley-Nordström invariant is given by the formula (d - sum(w)) / d.")
    print(f"The quasi-homogeneous degree is d = {d}.")
    print(f"The sum of the weights is sum(w) = {sum_of_weights}.")
    print(f"\nFinal Equation:")
    print(f"({d} - ({weights_str})) / {d} = {invariant}")

solve_crawley_nordstrom()