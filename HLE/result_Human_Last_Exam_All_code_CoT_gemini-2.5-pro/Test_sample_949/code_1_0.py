def solve():
    """
    Solves the parenthesis string analysis problem.

    The problem asks to evaluate six statements comparing functions of length L(x) and depth D(x)
    for all matching parenthesis pairs x in a string.

    1. sum(log L(x)) = O(sum(log D(x))): False.
       Consider a family of strings corresponding to a complete b-ary tree of depth d.
       Let b be fixed and d -> infinity.
       The LHS is asymptotically proportional to b^d * log(b), while the RHS is proportional to b^d.
       The ratio contains a log(b) factor, which can be made arbitrarily large by choosing a large b.
       Thus, no single constant C can satisfy the O-notation for all strings.

    2. sum(loglog L(x)) = O(sum(loglog D(x))): False.
       Consider a complete b-ary tree of fixed depth d=2 and let b -> infinity. This corresponds to
       a string like `( () () ... () )` with b inner pairs.
       The dominant term in the LHS sum involves log(log(b)), while the dominant term in the RHS is constant.
       The ratio is unbounded as b grows.

    3. sum(log^5 L(x)) = O(sum(log^5 D(x))): False.
       Similar to statement 2, the ratio of sums for the S_{d,b} family with d=2, b->inf
       grows like (log b)^5.

    4. sum(2^sqrt(log L(x))) = O(sum(2^sqrt(log D(x)))): False.
       Again, for the S_{d,b} family with d=2, b->inf, the ratio of sums
       grows like 2^(sqrt(log b)).

    5. sum(L(x)^0.1) = O(sum(D(x)^0.11)): True.
       The exponent on the right (0.11) is larger than on the left (0.1).
       When L(x) is much larger than D(x), it's because the corresponding tree is "wide,"
       meaning it has many nodes. The sum on the RHS accumulates over all these nodes,
       and the larger exponent provides enough power to bound the LHS. Attempts to
       construct counterexamples with various tree structures fail.

    6. sum(L(x)^(1/4)) = O(sum(D(x)^(1/2))): True.
       Similar to statement 5, the exponent on the right (0.5) is larger than on the left (0.25),
       making the statement true.

    The final result is a string of T/F values.
    """
    
    # The derived answers for statements 1 through 6.
    answers = ['F', 'F', 'F', 'F', 'T', 'T']
    
    result_string = "".join(answers)
    
    print(result_string)

solve()