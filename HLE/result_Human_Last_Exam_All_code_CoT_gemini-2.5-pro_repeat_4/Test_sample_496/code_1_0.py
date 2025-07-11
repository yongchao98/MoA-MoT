def solve():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4) \ X
    in degrees up to 100.
    """
    max_degree = 100

    # The Poincaré series for the equivariant cohomology ring A is (t^3 + t^6) / (1-t^4)^2.
    # We need to compute the sum of coefficients of this series for degrees up to 100.
    # The rank of H^j(BSO(4)) is m+1 if j=4m, and 0 otherwise.

    # Part 1: Contribution from the t^3 term.
    # This is the sum of ranks of H^j(BSO(4)) for j <= 100 - 3 = 97.
    # j must be a multiple of 4, so j=4m.
    # 4m <= 97 implies m <= 24.
    total_rank_1 = 0
    degree_shift_1 = 3
    max_m_1 = (max_degree - degree_shift_1) // 4
    for m in range(max_m_1 + 1):
        # The rank of the component of H*(BSO(4)) of degree 4m is m+1
        rank = m + 1
        total_rank_1 += rank

    # Part 2: Contribution from the t^6 term.
    # This is the sum of ranks of H^j(BSO(4)) for j <= 100 - 6 = 94.
    # j must be a multiple of 4, so j=4m.
    # 4m <= 94 implies m <= 23.
    total_rank_2 = 0
    degree_shift_2 = 6
    max_m_2 = (max_degree - degree_shift_2) // 4
    for m in range(max_m_2 + 1):
        # The rank of the component of H*(BSO(4)) of degree 4m is m+1
        rank = m + 1
        total_rank_2 += rank

    total_rank = total_rank_1 + total_rank_2

    print("The total rank is the sum of contributions from the two generators of the kernel of the restriction map.")
    print(f"The contribution from the degree 3 generator is the sum of ranks of H^*(BSO(4)) up to degree {max_degree - degree_shift_1}, which is: {total_rank_1}")
    print(f"The contribution from the degree 6 generator is the sum of ranks of H^*(BSO(4)) up to degree {max_degree - degree_shift_2}, which is: {total_rank_2}")
    print(f"The total rank of A in degrees * <= 100 is the final sum:")
    print(f"{total_rank_1} + {total_rank_2} = {total_rank}")

solve()