def solve_involutions():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # A is the number of involutions g in SU(4,q) where g is not I or -I.
    # Formula: A = q^4 * (q^2 - q + 1) * (q^2 + 1)
    q2 = q * q
    A = (q**4) * (q2 - q + 1) * (q2 + 1)

    # B is the number of elements g in SU(4,q) such that g^2 = -I.
    # Formula: B = q^5 * (q^2 - 1)^2 * (q^2 + 1) * (q^2 - q + 1)
    B = (q**5) * ((q2 - 1)**2) * (q2 + 1) * (q2 - q + 1)

    # The total number of involutions in PSU(4, q) is (A + B) / 2
    num_involutions = (A + B) // 2

    # The problem asks to output the numbers in the final equation.
    print(f"The number of involutions is determined by the formula (A + B) / 2.")
    print(f"For q = 997:")
    print(f"Number of 'true' SU(4,q) involutions (A) = {A}")
    print(f"Number of 'skew' SU(4,q) involutions (B) = {B}")
    print(f"Total number of involutions in PSU(4,997) = ({A} + {B}) / 2 = {num_involutions}")

solve_involutions()