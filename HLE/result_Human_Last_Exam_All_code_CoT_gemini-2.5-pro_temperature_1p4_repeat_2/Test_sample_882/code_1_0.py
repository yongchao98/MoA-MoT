import math

def solve_for_p():
    """
    Solves the game theory equilibrium equation E_D(p) = E_S(p) numerically.

    Returns:
        float: The equilibrium probability p.
    """
    N = 8
    M = 3 * N - 1
    M_plus_1 = 3 * N

    memo_h = {}

    def h_func(p, j):
        """
        Calculates a core component of the E_S(p) payoff function, which is an
        expectation involving a binomial distribution.
        H(p, j) = E_s[(1-j/N)^(M-s) * 1/(s+1)] where s ~ Bin(M, 1-p).
        """
        if (p, j) in memo_h:
            return memo_h[(p, j)]

        # The closed-form analytical expression for H(p, j) can be unstable when p is
        # very close to 1. A small tolerance is used to fall back to a direct,
        # more stable summation method in that edge case.
        if abs(p - 1.0) < 1e-9:
            # Fallback to direct summation for p very close to 1
            # When p=1, s=0, so this gives (1-j/N)^M
            total_sum = 0
            # s is the number of players choosing strategy S from M others.
            # s follows Binomial(M, 1-p).
            for s in range(M + 1):
                try:
                    # Compute in log-space to handle small numbers
                    log_term = (math.log(math.comb(M, s)) +
                                s * math.log(1 - p) +
                                (M - s) * math.log(p) +
                                (M - s) * math.log(1 - j / N) -
                                math.log(s + 1))
                    if log_term > -700:  # Avoid underflow on exp
                        total_sum += math.exp(log_term)
                except ValueError:
                    # Handles log(0) cases if p=0 or p=1 exactly
                    continue
            res = total_sum
        else:
            # Analytical closed-form for H(p, j)
            v = p * (1.0 - j / N)
            v_plus_1_minus_p = 1.0 - p * j / N
            num = v_plus_1_minus_p**M_plus_1 - v**M_plus_1
            den = M_plus_1 * (1.0 - p)
            res = num / den

        memo_h[(p, j)] = res
        return res

    def e_s(p):
        """
        Calculates the expected payoff for the Spreading Strategy (S_N).
        This is based on an inclusion-exclusion argument over the N races.
        """
        total = 0
        for j in range(1, N + 1):
            term = math.comb(N, j) * h_func(p, j)
            if j % 2 == 1:
                total += term
            else:
                total -= term
        return total

    def e_d(p):
        """
        Calculates the expected payoff for the Discrete Strategy (D).
        """
        if abs(p) < 1e-12:
            return 1.0

        num = 1.0 - (1.0 - p / N)**M_plus_1
        den = M_plus_1 * (p / N)
        return num / den

    def f(p):
        """
        The function whose root we want to find: f(p) = E_D(p) - E_S(p).
        """
        # Clear memoization for each new p value during search
        memo_h.clear()
        return e_d(p) - e_s(p)

    # Use bisection method to find the root p where f(p) = 0.
    low = 0.5
    high = 1.0

    # We know from the problem that f(1) < 0. For p near 0, f(p) > 0.
    # A root must exist in (0, 1).
    for _ in range(100):  # 100 iterations provides high precision
        mid = (low + high) / 2
        if mid == low or mid == high: # Precision limit reached
            break

        f_mid = f(mid)

        if f_mid < 0:
            high = mid
        else:
            low = mid

    return (low + high) / 2

# --- Main execution ---
p_solution = solve_for_p()
one_minus_p = 1 - p_solution
value_to_floor = 10000 * one_minus_p
result = math.floor(value_to_floor)

print(f"The equilibrium probability 'p' is found to be: {p_solution}")
print(f"The final equation is: floor(10000 * (1 - p))")
print(f"Substituting the value of p into the equation:")
print(f"floor(10000 * (1 - {p_solution:.6f}))")
print(f"= floor(10000 * {one_minus_p:.6f})")
print(f"= floor({value_to_floor:.6f})")
print(f"= {result}")