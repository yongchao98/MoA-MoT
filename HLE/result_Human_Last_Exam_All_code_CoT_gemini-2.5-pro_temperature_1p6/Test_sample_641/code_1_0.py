def solve_involutions():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    # In PSU(n, q), the number of involutions, for n=4, q=997, is given by
    # a specific formula derived from the group's structure.
    # The conditions for this formula are d=gcd(n, q+1)=2 and q=1 (mod 4),
    # which are met here (d=gcd(4,998)=2 and 997=1 (mod 4)).
    # The formula is: 1 + q^4 * (q^2 - q + 1) * (q^2 + 1).

    q = 997

    # Calculate the components of the formula
    q_pow_2 = q**2
    q_pow_4 = q**4

    term1 = q_pow_2 - q + 1
    term2 = q_pow_2 + 1

    # Calculate the final number of involutions
    num_involutions = 1 + q_pow_4 * term1 * term2

    # As requested, print the equation showing each evaluated number
    print(f"The number of involutions is calculated by the formula: 1 + q^4 * (q^2 - q + 1) * (q^2 + 1)")
    print(f"For q = 997, the equation with evaluated numbers is:")
    print(f"1 + {q_pow_4} * {term1} * {term2} = {num_involutions}")

solve_involutions()