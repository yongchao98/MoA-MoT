def solve():
    """
    Calculates the number of involutions in the group PSU(4, 997).
    """
    # Parameter for the group PSU(4, q)
    q = 997

    # The number of involutions in PSU(4, q) for q = 1 (mod 4) is given by the formula:
    # N = q^4 * (q^2 - q + 1) * (q^2 + 1) + 1
    # This formula is applicable as 997 % 4 == 1.
    
    # We use integer arithmetic to handle the large numbers involved.
    try:
        # Calculate each part of the formula
        q_4 = q**4
        term1 = q**2 - q + 1
        term2 = q**2 + 1
        one = 1

        # Compute the final result
        num_involutions = q_4 * term1 * term2 + one

        # Print the final equation with all numbers shown explicitly as requested.
        print(f"The number of involutions is calculated from the expression: q^4 * (q^2 - q + 1) * (q^2 + 1) + 1")
        print(f"For q = {q}, the equation is:")
        print(f"{q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1) + 1 = {num_involutions}")

    except Exception as e:
        print(f"An error occurred during calculation: {e}")

solve()