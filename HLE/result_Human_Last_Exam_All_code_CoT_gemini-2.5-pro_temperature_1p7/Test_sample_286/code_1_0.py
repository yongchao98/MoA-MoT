import math

def solve():
    """
    Solves the problem to find the smallest positive integer M.

    Let L = 1 - 10^{-M}. The interval for x_i is [-L, L].
    The condition is that for any sequence a_i in [0,1], there's a sequence x_i in [-L,L]
    such that |x_{i-1} - x_i| = a_i.

    This can be reformulated as a problem in discrepancy theory. Let x_k = x_0 + p_k,
    where p_k is a sum of signed a_i's. For the sequence {x_k} to fit into [-L, L],
    its span must be at most 2L. The span of {x_k} is the span of {p_k}.
    span({p_k}) = max_{j,k} |p_k - p_j|.

    So, for any {a_i}, we must be able to choose signs {s_i} such that this span is <= 2L.
    This leads to the condition: max_{a} (min_{s} (span({p_k}))) <= 2L.

    Let V_N be this max-min value for N steps. A known result from discrepancy theory is
    that for large N, V_N is approximately 2 - 2/N.
    We are given N = 100000.

    So, we need V_N <= 2 * (1 - 10^{-M}).
    2 - 2/N <= 2 - 2 * 10^{-M}
    -2/N <= -2 * 10^{-M}
    1/N >= 10^{-M}
    1/100000 >= 10^{-M}
    10^{-5} >= 10^{-M}
    -5 >= -M
    M >= 5

    The smallest positive integer M is 5.
    """
    N = 100000
    
    # We derived the inequality M >= log10(N)
    log_N = math.log10(N)
    
    # M must be the smallest positive integer satisfying M >= log10(N)
    M = math.ceil(log_N)

    # Output the explanation
    print(f"Let N = {N}.")
    print("The condition to be satisfied can be expressed using concepts from discrepancy theory.")
    print(f"Let V_N be the maximum possible minimal walk span for N steps.")
    print(f"For large N, V_N is approximately 2 - 2/N.")
    print(f"The condition is V_N <= 2 * (1 - 10**(-M)).")
    print(f"Substituting the approximation for V_N gives: 2 - 2/N <= 2 - 2 * 10**(-M).")
    print(f"This simplifies to 1/N >= 10**(-M).")
    print(f"Taking log10 of both sides: log10(1/N) >= -M, which is -log10(N) >= -M.")
    print(f"This is equivalent to M >= log10(N).")
    print(f"For N = {N}, log10(N) is {int(log_N)}.")
    print(f"So we need M >= {int(log_N)}.")
    print(f"The smallest positive integer M that satisfies this is {M}.")
    print(f"The final equation to solve for M is: 1 / {N} >= 10**(-M)")

solve()