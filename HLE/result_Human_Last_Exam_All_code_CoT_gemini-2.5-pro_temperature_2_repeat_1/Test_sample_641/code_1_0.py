def count_involutions_in_psu_4_q():
    """
    Calculates the number of involutions in PSU(4, q) for a given q.
    The number is given by the formula N = q^4 * (q^2 - q + 1) * (q^2 + 1).
    This script computes the value for q = 997.
    """
    q = 997

    # Calculate the terms of the formula
    q_squared = q**2
    q_power_4 = q**4
    term1 = q_squared - q + 1
    term2 = q_squared + 1

    # Calculate the final result
    result = q_power_4 * term1 * term2

    # Print the equation and the result, showing each number
    print(f"The number of involutions in PSU(4, {q}) is calculated as follows:")
    print(f"Formula: N = q^4 * (q^2 - q + 1) * (q^2 + 1)")
    print(f"Substituting q = {q}:")
    print(f"N = {q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1)")
    print(f"N = {q_power_4} * ({q_squared} - {q} + 1) * ({q_squared} + 1)")
    print(f"N = {q_power_4} * {term1} * {term2}")
    print(f"N = {result}")

# Run the calculation
count_involutions_in_psu_4_q()
