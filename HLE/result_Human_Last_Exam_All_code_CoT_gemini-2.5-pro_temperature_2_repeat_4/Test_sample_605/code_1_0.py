import sys

def calculate_crawley_nordstrom_invariant():
    """
    Calculates the Crawley-Nordström invariant for a Calabi-Yau link
    defined in a weighted projective space.
    """
    # 1. Define the weights for the weighted projective space.
    weights = [22, 29, 49, 50, 75]

    # 2. For a Calabi-Yau hypersurface, the degree of the defining
    #    polynomial 'd' must equal the sum of the weights 'Σw_i'.
    sum_of_weights = sum(weights)
    degree = sum_of_weights

    # 3. The Crawley-Nordström invariant is defined by the formula: (d - Σw_i) / d.
    #    This value is always zero for a Calabi-Yau hypersurface by definition.
    #    Let's calculate it step-by-step.
    
    # Check for division by zero, though unlikely in this context.
    if degree == 0:
        print("Error: The degree 'd' is zero, cannot divide by zero.", file=sys.stderr)
        return

    numerator = degree - sum_of_weights
    invariant = numerator / degree

    # 4. Print the steps of the calculation.
    print(f"The weights of the ambient space are: {weights}")
    print(f"The sum of the weights is: Σw_i = {' + '.join(map(str, weights))} = {sum_of_weights}")
    print("For a Calabi-Yau Link, the degree 'd' must equal the sum of the weights.")
    print(f"Therefore, the degree d = {degree}.")

    print("\nThe Crawley-Nordström invariant is calculated using the formula: (d - Σw_i) / d")
    print("Substituting the values into the formula gives the final equation:")
    print(f"({degree} - {sum_of_weights}) / {degree} = {invariant}")
    print(f"\nTherefore, the Crawley-Nordström invariant for the given Calabi-Yau Link is {invariant}.")

if __name__ == '__main__':
    calculate_crawley_nordstrom_invariant()