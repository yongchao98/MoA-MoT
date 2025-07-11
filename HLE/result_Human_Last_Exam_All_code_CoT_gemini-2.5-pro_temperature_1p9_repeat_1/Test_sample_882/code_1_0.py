import math

def solve():
    """
    Solves for the equilibrium probability p and computes the final value.
    """
    N = 8
    K = 3 * N  # Total number of players

    def g(p):
        """
        The function g(p) = pi_D(p) - pi_S(p). We need to find the root of this function.
        We handle the edge cases p=0 and p=1 using their limits to avoid division by zero.
        """
        if p == 0:
            # lim p->0 pi_D(p) = 1
            # lim p->0 pi_S(p) = (N/K) * ((N-1)/N)**K
            return 1.0 - (N / K) * math.pow((N - 1) / N, K)
        if p == 1:
            # lim p->1 pi_D(p) = (N/K) * (1 - ((N-1)/N)**K)
            # lim p->1 pi_S(p) = N - 1
            return (N / K) * (1.0 - math.pow((N - 1) / N, K)) - (N - 1)

        # Payoff for Discrete Strategy
        term_d = math.pow(1.0 - p / N, K)
        pi_D = (N / (p * K)) * (1.0 - term_d)

        # Payoff for Splitting Strategy
        term_s1 = math.pow((N - 1 + p) / N, K)
        term_s2 = math.pow(p, K)
        pi_S = (N / ((1 - p) * K)) * (term_s1 - term_s2)
        
        return pi_D - pi_S

    # Bisection method to find the root of g(p)
    low = 0.0
    high = 1.0
    
    # Check if a root exists in the interval
    if g(low) * g(high) >= 0:
        print("Bisection method may fail: g(0) and g(1) do not have opposite signs.")
        return

    for _ in range(100):  # 100 iterations for high precision
        mid = (low + high) / 2.0
        if g(mid) == 0:
            p = mid
            break
        elif g(low) * g(mid) < 0:
            high = mid
        else:
            low = mid
    
    p = (low + high) / 2.0
    
    # The final equation with the calculated value of p
    # To avoid clutter, we show the main structure without the full expansion of powers
    equation_str = (f"Equation with p = {p:.6f}:\n"
                    f"({N}/({K}*{p:.6f})) * [1 - (1 - {p:.6f}/{N})^{K}] = "
                    f"({N}/({K}*(1 - {p:.6f}))) * [(( {N - 1} + {p:.6f})/{N})^{K} - {p:.6f}^{K}]")

    print(equation_str)
    
    # Calculate the final value
    result = math.floor(10000 * (1.0 - p))
    print("\nThe value of p is approximately {:.6f}".format(p))
    print("The result of floor(10000 * (1-p)) is:")
    print(result)
    
solve()
<<<2991>>>