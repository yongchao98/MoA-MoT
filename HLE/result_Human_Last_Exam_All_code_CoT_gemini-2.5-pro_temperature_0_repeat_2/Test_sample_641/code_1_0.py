def solve():
    """
    Calculates the number of involutions in PSU(4, 997).
    """
    q = 997

    # The number of involutions in PSU(4, q) for d=gcd(4,q+1)=2 is given by:
    # (Number of non-central involutions in SU(4,q) + Number of square roots of -I in SU(4,q)) / 2
    # Let's calculate these two components.

    # Component A: Number of non-central involutions in SU(4, q).
    # These are involutions with a 2-dimensional (-1)-eigenspace.
    # The number is given by the formula N(4,2) = q^4 * (q^2 - q + 1) * (q^2 + 1).
    q2 = q * q
    q4 = q2 * q2
    term_factor1 = q2 - q + 1
    term_factor2 = q2 + 1
    num_non_central_involutions_su = q4 * term_factor1 * term_factor2

    # Component B: Number of square roots of -I in SU(4, q).
    # For q=997 (which is 1 mod 4), this is |U(4,q)| / |GL(2,q)|.
    # This simplifies to q^5 * (q+1)^3 * (q^2 - q + 1) * (q^2 + 1).
    q5 = q4 * q
    q_plus_1 = q + 1
    q_plus_1_cubed = q_plus_1 * q_plus_1 * q_plus_1
    num_sqrt_neg_one_su = q5 * q_plus_1_cubed * term_factor1 * term_factor2

    # The total number of involutions in PSU(4, q) is the sum of these two quantities, divided by 2.
    total_involutions = (num_non_central_involutions_su + num_sqrt_neg_one_su) // 2

    print(f"The number of involutions in PSU(4, {q}) is calculated as (A + B) / 2, where:")
    print(f"A = Number of non-central involutions in SU(4, {q}) = {num_non_central_involutions_su}")
    print(f"B = Number of square roots of -I in SU(4, {q}) = {num_sqrt_neg_one_su}")
    print(f"The final calculation is:")
    print(f"({num_non_central_involutions_su} + {num_sqrt_neg_one_su}) / 2 = {total_involutions}")

solve()