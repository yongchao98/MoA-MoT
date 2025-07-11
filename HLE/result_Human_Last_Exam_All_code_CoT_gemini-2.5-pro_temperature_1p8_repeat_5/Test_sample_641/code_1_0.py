def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    # The field parameter q
    q = 997

    # The formula for the number of involutions in PSU(4, q) for q=1 (mod 4) is:
    # N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2

    # Calculate the components of the formula
    q_pow_2 = q**2
    q_pow_4 = q_pow_2**2
    term1 = q_pow_2 - q + 1
    term2 = q_pow_2 + 1

    # The number of non-central involutions in SU(4, q)
    num_non_central_involutions_su = q_pow_4 * term1 * term2

    # The number of involutions in PSU(4, q)
    num_involutions_psu = num_non_central_involutions_su // 2

    # Print the equation with the calculated numbers and the final result
    print("The number of involutions in PSU(4, q) is given by the formula:")
    print("N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2")
    print(f"For q = {q}:")
    print(f"N = ({q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1)) / 2")
    print(f"N = ({q_pow_4} * {term1} * {term2}) / 2")
    print(f"N = {num_non_central_involutions_su} / 2")
    print(f"N = {num_involutions_psu}")

solve()