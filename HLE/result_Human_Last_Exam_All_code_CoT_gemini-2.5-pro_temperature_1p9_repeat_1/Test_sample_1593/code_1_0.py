import math

def solve():
    """
    Calculates the query complexity for the two specified regimes.
    
    The analysis shows that the overall complexity is governed by the minimum of
    two sorting strategies: a direct comparison sort with complexity O(N log N)
    and a radix sort with complexity O(N * L / log N).
    
    The final complexity is therefore O(N * min(log N, L / log N)).
    """
    
    # We represent the complexity Theta(N * F(N)) as (a, b, c)
    # Theta(sqrt(N^a * (log N)^b * (log log N)^c))
    # Squaring gives N^2 * F(N)^2 = N^a * (log N)^b * (log log N)^c
    
    results = []

    # Regime 1: N = 2^sqrt(L)  => L = (log N)^2
    # We compare log(N) with L/log(N) = (log N)^2 / log N = log N
    # The two terms are equal.
    # So, the complexity is Theta(N * log N).
    # N^2 * (log N)^2 = N^a * (log N)^b * (log log N)^c
    a1 = 2
    b1 = 2
    c1 = 0
    regime1_complexity = (a1, b1, c1)
    results.append(regime1_complexity)

    # Regime 2: N = 2^((log L)^2) => L = 2^sqrt(log N)
    # We compare log(N) with L/log(N) = 2^sqrt(log N) / log N
    # Let x = log N. We compare x with 2^sqrt(x)/x.
    # This is equivalent to comparing x^2 with 2^sqrt(x).
    # Taking log2 of both, we compare 2*log2(x) with sqrt(x).
    # Asymptotically, sqrt(x) grows faster than 2*log2(x).
    # Therefore, 2^sqrt(x) grows faster than x^2 for large x.
    # So, L/log(N) is larger than log(N) asymptotically.
    # The minimum term is log(N).
    # So, the complexity is Theta(N * log N).
    # N^2 * (log N)^2 = N^a * (log N)^b * (log log N)^c
    a2 = 2
    b2 = 2
    c2 = 0
    regime2_complexity = (a2, b2, c2)
    results.append(regime2_complexity)

    # Format the output string
    output_str = f"({results[0][0]},{results[0][1]},{results[0][2]}),({results[1][0]},{results[1][1]},{results[1][2]})"
    print(output_str)

solve()