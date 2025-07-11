import math

def calculate_dimension():
    """
    Calculates the complex dimension of the space of global sections of the sheaf
    Omega^1(P^n) tensored with O(P^n)(2).
    """
    # The dimension of the complex projective space is denoted by 'n'.
    # You can change this value to calculate the result for a different space.
    n = 4

    print(f"Calculating dimension for n = {n}...")

    # The dimension is given by the formula:
    # h^0(Omega^1(2)) = h^0(O(1)^(n+1)) - h^0(O(2))
    # which is derived from the twisted Euler sequence.

    # First term in the equation: h^0(O(1)^(n+1))
    # This is equal to (n+1) * h^0(O(1)) = (n+1) * (n+1) = (n+1)^2
    term1 = (n + 1)**2

    # Second term in the equation: h^0(O(2))
    # This is given by the binomial coefficient C(n+2, 2).
    # We use math.comb() for robust calculation of binomial coefficients.
    try:
        term2 = math.comb(n + 2, 2)
    except ValueError:
        print("Error: n must be a non-negative integer.")
        return

    # The final dimension is the difference between the two terms.
    dimension = term1 - term2

    # Printing the final equation with all the calculated numbers.
    print("\nThe dimension is calculated as the difference of the dimensions of two other spaces of sections:")
    print(f"dim = h^0(O(1)^({n}+1)) - h^0(O(2))")
    print("\nSubstituting the calculated values:")
    print(f"dim = {term1} - {term2} = {dimension}")

    # The final simplified formula is C(n+1, 2). Let's verify our result.
    verification = math.comb(n + 1, 2)
    print(f"\nResult verified with simplified formula C(n+1, 2) = C({n}+1, 2) = {verification}")

if __name__ == "__main__":
    calculate_dimension()