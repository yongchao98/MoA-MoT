import math

def solve():
    """
    Solves for the mixed strategy equilibrium probability p and computes the final value.
    """
    # Step 1: Set up constants for N=8
    N = 8
    M = 3 * N - 1  # Total other players
    q_N = 1.0 - 1.0 / N # (1 - 1/N)

    # Step 2: Define the function f(p) = U_discrete(p) - U_spread(p) whose root we want to find.
    def f(p):
        """
        Calculates the difference in expected payoffs between playing discrete and spread strategies.
        """
        if p < 0 or p > 1:
            return float('nan')
        if p == 0:
            # U_discrete(0)=1, U_spread(0)=1/3.
            return 1.0 - 1.0/3.0
        if p == 1:
            # U_discrete(1) = (1 - (1-1/N)^(M+1)) / ((M+1)/N)
            # U_spread(1) = N * (1-1/N)^M
            ud_1 = (1.0 - q_N**(M + 1)) / ((M + 1.0) / N)
            us_1 = N * (q_N**M) / (M - M + 1.0)
            return ud_1 - us_1

        expected_ud = 0.0
        expected_us = 0.0
        
        # Calculate the expectation over K ~ Binom(M, p)
        for K in range(M + 1):
            # Binomial probability P(K), computed using logs for numerical stability
            log_p_k = (math.lgamma(M + 1) - math.lgamma(K + 1) - math.lgamma(M - K + 1) +
                       K * math.log(p) + (M - K) * math.log(1 - p))
            
            prob_k = math.exp(log_p_k)

            # Payoff for discrete strategy, given K
            ud_k = (1.0 - q_N**(K + 1)) / ((K + 1.0) / N)
            
            # Payoff for spread strategy, given K
            us_k = N * (q_N**K) / (M - K + 1.0)
            
            expected_ud += prob_k * ud_k
            expected_us += prob_k * us_k
            
        return expected_ud - expected_us

    # Step 3: Use bisection method to find the root p of f(p)=0
    low_p = 0.0
    high_p = 1.0

    # f(0) is positive, f(1) is negative, so a root must exist in (0,1).
    for _ in range(100):  # 100 iterations for high precision
        mid_p = (low_p + high_p) / 2.0
        if mid_p == low_p or mid_p == high_p: # Precision limit reached
            break
        
        f_mid = f(mid_p)
        
        if f_mid > 0:
            low_p = mid_p
        else:
            high_p = mid_p

    p_star = (low_p + high_p) / 2.0

    # Step 4: Calculate the final value as requested and print the breakdown.
    one_minus_p = 1.0 - p_star
    value = 10000.0 * one_minus_p
    final_answer = math.floor(value)

    print("The final requested quantity is floor(10000 * (1 - p))")
    print(f"For N=8, the calculated equilibrium probability p is: {p_star:.8f}")
    print(f"The number 10000 is a constant from the problem.")
    print(f"The number 1 is a constant from the problem.")
    print(f"The equation with the value of p substituted is: floor(10000 * (1 - {p_star:.8f}))")
    print(f"This simplifies to: floor(10000 * {one_minus_p:.8f})")
    print(f"The intermediate value is: floor({value:.8f})")
    print(f"The final integer result is: {final_answer}")

solve()