def solve():
    """
    Calculates the total rank of the equivariant cohomology ring A
    in degrees up to 100.
    """
    max_degree = 100

    # The ranks are non-zero only for degrees i = 4k and i = 4k+3.
    # The rank is k+1 in both cases.

    # Sum ranks for degrees i = 4k <= 100
    # 4k <= 100  =>  k <= 25
    # The loop for k goes from 0 to 25.
    k_max_4k = max_degree // 4
    sum_ranks_4k = 0
    for k in range(k_max_4k + 1):
        rank = k + 1
        sum_ranks_4k += rank

    # Sum ranks for degrees i = 4k+3 <= 100
    # 4k <= 97  =>  k <= 24.25  => k <= 24
    # The loop for k goes from 0 to 24.
    k_max_4k_plus_3 = (max_degree - 3) // 4
    sum_ranks_4k_plus_3 = 0
    for k in range(k_max_4k_plus_3 + 1):
        rank = k + 1
        sum_ranks_4k_plus_3 += rank

    total_rank = sum_ranks_4k + sum_ranks_4k_plus_3
    
    print(f"The total rank is the sum of ranks from two series of terms.")
    print(f"Sum of ranks for degrees of the form 4k: {sum_ranks_4k}")
    print(f"Sum of ranks for degrees of the form 4k+3: {sum_ranks_4k_plus_3}")
    print(f"Total rank up to degree {max_degree}: {sum_ranks_4k} + {sum_ranks_4k_plus_3} = {total_rank}")

solve()