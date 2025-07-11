import math

def solve_for_M():
    """
    Solves the mathematical problem to find the smallest integer M.

    The problem requires finding the smallest positive integer M such that for any
    a_i in [0,1], we can find x_i in [-S, S] where S = 1 - 10**(-M),
    satisfying |x_{i-1} - x_i| = a_i.

    This can be rephrased as: for any a_i, there must exist a choice of signs
    eps_i in {-1, 1} and a starting point x_0 in [-S, S] such that all
    x_k = x_0 + sum_{i=1 to k} (eps_i * a_i) are in [-S, S].

    This is possible if and only if the range of the partial sums,
    R = max_k(y_k) - min_k(y_k) where y_k = sum_{i=1 to k} (eps_i * a_i),
    is at most 2*S.

    So, for any sequence {a_i}, we need to be able to choose {eps_i} such that
    the range R <= 2*S. This means that the "worst-case minimal range" must
    satisfy the inequality.
    
    Let K = sup_{a_i} min_{eps_i} (Range). The condition is K <= 2*S.
    
    A key result from discrepancy theory states that K = 1. This worst-case value
    is achieved for the sequence a_i = 1 for all i. For this sequence, a greedy
    choice of signs yields a range of 1.

    Therefore, the problem reduces to solving the inequality for M:
    1 <= 2 * (1 - 10**(-M))
    """

    print("Step 1: State the core inequality derived from the problem's conditions.")
    print("The problem requires that the maximum possible value of the minimum range of partial sums is less than or equal to the effective width of the interval for the sequence x_i.")
    print("Let K = sup_{a_i} min_{eps_i} (range of partial sums).")
    print("Let S = 1 - 10**(-M). The condition is K <= 2 * S.")
    
    print("\nStep 2: Determine the value of K.")
    print("Based on results in discrepancy theory, the value of K is 1. This corresponds to the case where a_i = 1 for all i.")
    
    print("\nStep 3: Formulate the final inequality with K=1.")
    equation_str_1 = "1 <= 2 * (1 - 10**(-M))"
    print(f"The inequality to solve is: {equation_str_1}")

    print("\nStep 4: Solve the inequality for M.")
    print("Divide by 2:")
    print("1 / 2 <= 1 - 10**(-M)")
    print("0.5 <= 1 - 10**(-M)")
    print("\nRearrange the terms to isolate 10**(-M):")
    print("10**(-M) <= 1 - 0.5")
    equation_str_2 = "10**(-M) <= 0.5"
    print(f"{equation_str_2}")
    
    print("\nTake the base-10 logarithm of both sides:")
    print("log10(10**(-M)) <= log10(0.5)")
    print("-M <= log10(0.5)")
    print("-M <= -log10(2)")
    
    print("\nMultiply by -1 (and flip the inequality sign):")
    equation_str_3 = "M >= log10(2)"
    print(f"{equation_str_3}")

    log10_of_2 = math.log10(2)
    print(f"\nStep 5: Calculate the numerical value and find the smallest integer M.")
    print(f"The value of log10(2) is approximately {log10_of_2:.6f}.")
    
    # M must be the smallest positive integer satisfying the condition.
    M = math.ceil(log10_of_2)
    print(f"Since M must be a positive integer, the smallest possible value for M is ceil({log10_of_2:.6f}), which is {M}.")
    
    print("\nFinal Answer:")
    print(f"M = {M}")

solve_for_M()