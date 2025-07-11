import math

def solve():
    """
    Solves the problem by explaining the condition and calculating the limit.
    """
    
    # Part 1: The condition for a pair (p,q) to be good.
    # p, q are coprime integers such that 1 < p,q <= n.
    # The pair (p,q) is "good" if the graph G with vertices {1,...,n} and
    # edges between i and j if |i-j| in {p,q} is connected.

    condition_explanation = """
A pair (p,q) is good if and only if the corresponding graph is connected.
This graph is disconnected if and only if, for p or q (let's assume p < q without loss of generality),
one of the vertex sets C_r = {x | 1<=x<=n, x = r (mod p)} is 'isolated' by q-swaps.
A set C_r is isolated if for every x in C_r, x+q > n and x-q < 1.

A simpler sufficient condition for (p,q) to be good is: p + q <= n.
All pairs that do not satisfy p+q <= n are not necessarily bad. The pairs that
are not good form a very small subset of all possible pairs, especially as n grows large.
    """
    print("--- Condition for a (p,q) pair to be good ---")
    print(condition_explanation)

    # Part 2: The limit of the probability Pr(n).
    # Pr(n) is the probability that selecting p,q randomly from {2,...,n},
    # the pair is coprime and good.
    # Total pairs = (n-1)^2
    # Number of coprime pairs is asymptotically (n-1)^2 * (6 / pi^2)
    # The number of "bad" pairs is asymptotically negligible, i.e., o(n^2).
    # This is because the conditions for a pair to be bad are very restrictive.
    # It can be shown that if p+q <= n, the pair is always good. Bad pairs must
    # have p+q > n. The number of such bad pairs is small compared to n^2.
    # So, the number of good coprime pairs is asymptotically equal to the
    # total number of coprime pairs.
    
    # Pr(n) = (Number of good coprime pairs) / (Total pairs)
    # lim Pr(n) = lim ( (6/pi^2)*(n-1)^2 / (n-1)^2 ) = 6 / pi^2
    
    numerator = 6
    # In Python, math.pi is an approximation of pi.
    denominator_symbolic = "pi^2"
    denominator_numeric = math.pi**2
    
    limit_value = numerator / denominator_numeric
    
    print("\n--- Calculation of the Limit ---")
    print("The limit of Pr(n) as n -> infinity is given by the equation:")
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"lim Pr(n) = {numerator} / {denominator_symbolic}")
    print(f"The numerical value of the limit is {numerator} / {denominator_numeric:.7f} ...")
    print(f"Limit value: {limit_value}")
    
    final_answer = limit_value
    return final_answer
    
if __name__ == '__main__':
    # The final answer is the exact value, not the code itself.
    # We execute the code to derive and print the final answer.
    solve()
    # The actual numerical value is 6 / (pi^2)
    final_limit = 6 / (math.pi**2)
    print(f"\n<<<6/pi^2>>>")
