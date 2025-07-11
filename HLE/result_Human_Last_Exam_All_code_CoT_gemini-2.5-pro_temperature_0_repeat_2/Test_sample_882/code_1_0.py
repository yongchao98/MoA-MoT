import math

def combinations(n, k):
    """Calculates the number of combinations nCk."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def main():
    """
    Solves the game theory problem to find the mixed strategy equilibrium probability.
    """
    N = 8
    num_players = 3 * N
    n_others = num_players - 1

    # 1. Calculate E_D(1) and E_A(1)
    # E_D(1) is the payoff for a player using the discrete strategy when all n_others
    # also use the discrete strategy.
    # Formula: (1 - (1 - 1/N)^(n_others + 1)) / ((n_others + 1) * (1/N))
    e_d_at_1 = (1 - (1 - 1/N)**(n_others + 1)) / ((n_others + 1) * (1/N))

    # E_A(1) is the payoff for a player using the alternative (splitting) strategy
    # when all n_others use the discrete strategy.
    # Formula: sum_{j=1 to N} (-1)^(j-1) * C(N,j) * (1 - j/N)^n_others
    e_a_at_1 = 0
    for j in range(1, N + 1):
        term = ((-1)**(j - 1)) * combinations(N, j) * ((N - j) / N)**n_others
        e_a_at_1 += term

    # 2. Calculate the derivatives E'_D(1) and E'_A(1)
    
    # To calculate E'_D(1), we use the formula n*(G(n) - G(n-1)), where n=n_others.
    def G(k, N_val):
        if k < 0: return 0
        # Payoff for discrete player against k discrete players
        return (N_val * (1 - (1 - 1/N_val)**(k + 1))) / (k + 1)

    e_prime_d_at_1 = n_others * (G(n_others, N) - G(n_others - 1, N))

    # To calculate E'_A(1), we use the formula for the derivative of E_A(p) at p=1.
    # Formula: -(n/N) * sum_{j=1 to N} (-1)^(j-1) * j * C(N,j) * (1-j/N)^(n-1)
    sum_term_for_e_prime_a = 0
    for j in range(1, N + 1):
        term = ((-1)**(j - 1)) * j * combinations(N, j) * ((N - j) / N)**(n_others - 1)
        sum_term_for_e_prime_a += term
    e_prime_a_at_1 = -(n_others / N) * sum_term_for_e_prime_a

    # 3. Calculate delta = 1 - p
    numerator = e_d_at_1 - e_a_at_1
    denominator = e_prime_d_at_1 - e_prime_a_at_1
    delta = numerator / denominator
    
    p = 1 - delta

    # 4. Output the results
    print("For N = 8, the equilibrium is found by solving E_D(p) = E_A(p).")
    print("We approximate p near 1 using a first-order Taylor expansion.")
    print("This gives the equation for delta = 1 - p:")
    print(f"delta = (E_D(1) - E_A(1)) / (E'_D(1) - E'_A(1))")
    print("\nCalculating the components:")
    print(f"E_D(1) = {e_d_at_1:.6f}")
    print(f"E_A(1) = {e_a_at_1:.6f}")
    print(f"E'_D(1) = {e_prime_d_at_1:.6f}")
    print(f"E'_A(1) = {e_prime_a_at_1:.6f}")
    
    print("\nPlugging the numbers into the equation:")
    print(f"1 - p = ({e_d_at_1:.6f} - {e_a_at_1:.6f}) / ({e_prime_d_at_1:.6f} - ({e_prime_a_at_1:.6f}))")
    print(f"1 - p = {numerator:.6f} / {denominator:.6f}")
    print(f"1 - p = {delta:.6f}")
    
    p_six_sig = f"{p:.6f}"
    print(f"\nThe probability p is approximately {p_six_sig}.")

    final_answer = math.floor(10000 * delta)
    print(f"\nThe final calculation is floor(10000 * (1-p)) = floor(10000 * {delta:.6f}) = {final_answer}")

if __name__ == "__main__":
    main()