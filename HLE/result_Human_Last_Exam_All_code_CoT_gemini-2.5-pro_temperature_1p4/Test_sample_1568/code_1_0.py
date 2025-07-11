def solve_infinite_product():
    """
    This function derives and prints the formula for the infinite product
    prod_{n=3 to inf} (1 - z^3/n^3).
    """

    # The problem asks for the product from n=3 to infinity.
    # We start with the known identity for the product from n=1.
    print("To solve this, we use a known product identity related to the Gamma function.")
    print("The infinite product from n=1 is given by:")
    print("  prod_{n=1 to inf} (1 - z^3/n^3) = 1 / [ Gamma(1-z) * Gamma(1-z*w) * Gamma(1-z*w^2) ]")
    print("  where w = exp(2*pi*i/3) is a complex cube root of unity, and Gamma is the Gamma function.\n")

    # The requested product starts from n=3, so we need to divide by the first two terms (n=1, n=2).
    n1 = 1
    n2 = 2
    
    # The cubes of n1 and n2 appear in the formula. Let's call them out as the "numbers" in the equation.
    term1_denominator = n1**3
    term2_denominator = n2**3

    print(f"The product we want starts from n=3. Therefore, we must divide the full product by the terms for n={n1} and n={n2}.")
    print(f"  - The term for n={n1} is: (1 - z^3 / {term1_denominator})")
    print(f"  - The term for n={n2} is: (1 - z^3 / {term2_denominator})\n")

    print("After dividing the full product by these two terms and simplifying, we get the final expression:")

    # Constructing the final formula as a string for printing.
    # The division results in the denominator of the n=2 term's factor (which is 8) moving to the numerator.
    final_numerator = term2_denominator
    final_denominator = f"({term1_denominator} - z^3) * ({term2_denominator} - z^3) * Gamma(1-z) * Gamma(1-z*w) * Gamma(1-z*w^2)"

    # Print the final formula in a visually clear format.
    print("\nResult:")
    print(f"          {final_numerator}")
    print("--------------------------------------------------------------------------")
    print(f"{final_denominator}\n")

solve_infinite_product()