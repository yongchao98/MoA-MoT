def solve():
    """
    Calculates the number of involutions in the projective special unitary group PSU(4, 997).
    """
    q = 997

    # According to the formulas for the number of pre-images of involutions in SU(4, q):
    # N1 is the number of non-central elements g such that g^2 = I.
    n1 = (q**4) * (q**2 + 1) * (q**2 - q + 1)

    # N2 is the number of elements g such that g^2 = -I.
    n2 = (q**2) * (q - 1) * (q**3 - 1)

    # The total number of involutions in PSU(4, q) is half the sum of these counts,
    # as each pair {g, -g} in SU(4,q) maps to a single involution in PSU(4,q).
    total_sum = n1 + n2
    total_involutions = total_sum // 2

    # Print the equation with the calculated numbers as requested.
    print(f"The number of involutions is calculated as (N1 + N2) / 2.")
    print(f"For q = {q}:")
    print(f"N1 = {q}^4 * ({q}^2 + 1) * ({q}^2 - {q} + 1)")
    print(f"   = {n1}")
    print(f"N2 = {q}^2 * ({q} - 1) * ({q}^3 - 1)")
    print(f"   = {n2}")
    print(f"\nFinal Equation:")
    print(f"Number of involutions = ({n1} + {n2}) / 2")
    print(f"                       = {total_sum} / 2")
    print(f"                       = {total_involutions}")

solve()