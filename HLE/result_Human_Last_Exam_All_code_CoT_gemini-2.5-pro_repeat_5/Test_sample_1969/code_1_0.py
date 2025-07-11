def print_shapley_formula(n):
    """
    For a given number of people n, this function calculates and prints the
    exact formula for the Shapley value c_k (the fair share for person p_k).
    The formula for c_k is a quadratic in k: A*k - B*k^2.
    """
    if not isinstance(n, int) or n <= 1:
        print("Error: Input 'n' must be an integer greater than 1.")
        return

    # S1 is the sum of the first n integers: 1 + 2 + ... + n
    S1 = n * (n + 1) // 2
    # S2 is the sum of the first n squared integers: 1^2 + 2^2 + ... + n^2
    S2 = n * (n + 1) * (2 * n + 1) // 6

    # The derived formula for the Shapley value c_k is:
    # c_k = k * S1 * (S2 + S1**2 - S1*k)
    # We can express this as a polynomial in k: c_k = A*k - B*k^2
    # where A = (S1*S2 + S1**3) and B = (S1**2)

    # Calculate the coefficients A and B
    coeff_A = S1 * S2 + S1**3
    coeff_B = S1**2

    print(f"For a group of n={n} people, the fair division c_k for person p_k is given by:")
    # The prompt asks to output each number in the final equation.
    # First, we show the structure with S1 and S2.
    print(f"\nc_k = k * {S1} * ({S2} + {S1}^2 - {S1} * k)")
    
    # Then, we show the simplified quadratic formula with calculated coefficients.
    print("\nThis simplifies to the quadratic formula:")
    print(f"c_k = {coeff_A} * k - {coeff_B} * k^2")

# --- User Execution ---
# You can execute this code. It will prompt you to enter a value for n.
try:
    n_value = int(input("Enter the total number of people, n (where n > 1): "))
    print_shapley_formula(n_value)
except ValueError:
    print("Invalid input. Please enter an integer.")
