import math

def calculate_torsion_rank(k_orig, n_orig):
    """
    Calculates the rank of the torsion subgroup of the integral cohomology
    of the real Grassmannian Gr(k, n).
    """

    # Use the duality Gr(k, n) = Gr(n-k, n) to simplify.
    # The formula is stated most easily for k <= n-k.
    k = min(k_orig, n_orig - k_orig)
    n = n_orig
    m = n - k # m = max(k_orig, n_orig - k_orig)

    print(f"The space is the Grassmannian Gr({k_orig}, {n_orig}).")
    print(f"Using the duality Gr({k_orig}, {n_orig}) ~ Gr({k}, {n}), we will use parameters k={k} and n={n}.")
    print("The rank of the torsion subgroup is calculated using Fomenko's formula:")
    print("Rank = Sum_{j=1 to k} N(n-k, j)")
    print("where N(m, j) is the count of integers s such that m < s <= m+j and C(s-1, j) is odd.")
    print(f"For this problem, k = {k} and m = n-k = {m}.\n")

    total_rank = 0
    N_values = []
    
    # j ranges from 1 to k
    for j in range(1, k + 1):
        N_j = 0
        print(f"--- Calculating N_{j} for j={j} ---")
        # s ranges from m+1 to m+j
        for s in range(m + 1, m + j + 1):
            # C(n, k) is 0 if k > n
            if s - 1 < j:
                comb = 0
            else:
                comb = math.comb(s - 1, j)

            # Check if the binomial coefficient is odd
            if comb % 2 == 1:
                print(f"s={s}: C({s-1}, {j}) = {comb}, which is odd. This contributes to N_{j}.")
                N_j += 1
            else:
                print(f"s={s}: C({s-1}, {j}) = {comb}, which is even.")
        
        print(f"The value for N_{j} is {N_j}.\n")
        N_values.append(N_j)
        total_rank += N_j

    # Build the final equation string
    rank_sum_symbols = " + ".join([f"N_{i+1}" for i in range(k)])
    rank_sum_values = " + ".join(map(str, N_values))

    print("--- Final Calculation ---")
    print(f"The total rank of the torsion subgroup is the sum of these values:")
    print(f"Rank = {rank_sum_symbols} = {rank_sum_values} = {total_rank}")


if __name__ == "__main__":
    # Parameters for the space of 3-subspaces of R^5
    k_param = 3
    n_param = 5
    calculate_torsion_rank(k_param, n_param)