def solve():
    """
    Solves the true/false puzzle about parenthesis strings.

    The reasoning for each statement is as follows:

    Let logp(z) = max(1, log(z)) to handle the well-definedness requirement.
    The primary test case is a "wide, shallow" string S_k = `(()()...())` with k inner pairs `()`
    and one outer pair. This string has k+1 total pairs. For this string, the sums are often
    dominated by the k inner pairs, which have simple L=2, D=1 values.

    1. sum logp(L) = O(sum logp(D)):
       For S_k, LHS is `logp(2k+2) + k*logp(2)` which is ~k.
       RHS is `logp(2) + k*logp(1)` which is ~k. Ratio is bounded.
       For other string structures, this appears to hold as well.
       Result: True.

    2. sum log(logp(L)) = O(sum log(logp(D))):
       For S_k, `logp(1)=1` and `logp(2)=1`. `log(1)=0`.
       LHS is `log(logp(2k+2)) + k*log(logp(2)) = log(log(2k+2)) + k*0`, which grows as log(log(k)).
       RHS is `log(logp(2)) + k*log(logp(1)) = 0 + k*0 = 0`.
       The ratio is unbounded.
       Result: False.

    3. sum logp(L)^5 = O(sum logp(D)^5):
       For S_k, LHS is `logp(2k+2)^5 + k*logp(2)^5 = (log(2k+2))^5 + k`, which is ~k.
       RHS is `logp(2)^5 + k*logp(1)^5 = 1 + k`, which is ~k.
       The ratio is bounded.
       Result: True.

    4. sum 2^sqrt(logp(L)) = sum 2^O(sqrt(logp(D))):
       This means `sum 2^sqrt(logp L) <= sum 2^(C*sqrt(logp D))` for some constant C.
       For S_k, LHS is `2^sqrt(log(2k+2)) + k*2^sqrt(1) = 2^sqrt(log(2k+2)) + 2k`.
       RHS is `(k+1) * 2^C`. We need to show `2k + small_term <= C'(k+1)`.
       This is true for a large enough C (e.g., C>=1, making C'>=2).
       Result: True.

    5. sum L^0.1 = O(sum D^0.11):
       For S_k, LHS is `(2k+2)^0.1 + k*2^0.1`, which is `~k * 2^0.1`.
       RHS is `2^0.11 + k*1^0.11`, which is `~k`.
       The ratio is `~2^0.1`, a constant. The larger exponent on D in the RHS helps ensure this holds for other structures too.
       Result: True.

    6. sum L^(1/4) = O(sum D^(1/2)):
       Similar to statement 5, the exponent on D (0.5) is larger than on L (0.25).
       For S_k, LHS is `(2k+2)^0.25 + k*2^0.25`, which is `~k * 2^0.25`.
       RHS is `2^0.5 + k*1^0.5`, which is `~k`.
       The ratio is `~2^0.25`, a constant.
       Result: True.
    """
    
    # The final deduced answer based on the analysis above.
    answer = "TFTTTT"
    print(answer)

solve()