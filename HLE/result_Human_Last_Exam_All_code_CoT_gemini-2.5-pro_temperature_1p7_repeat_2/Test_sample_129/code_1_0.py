import collections

def solve_boolean_expressions():
    """
    Calculates the number of valid boolean expressions of length 5
    using dynamic programming based on an unambiguous grammar.
    """
    max_len = 5
    # c: total expressions (Expr)
    # a: and-expressions (AndExpr)
    # m: not-expressions (NotExpr)
    # p: primary expressions (Primary)
    c = collections.defaultdict(int)
    a = collections.defaultdict(int)
    m = collections.defaultdict(int)
    p = collections.defaultdict(int)

    for n in range(1, max_len + 1):
        # Calculate p(n)
        p[n] = (2 if n == 1 else 0) + c.get(n - 2, 0)
        
        # Calculate m(n)
        m[n] = p[n] + m.get(n - 1, 0)
        
        # Calculate a(n)
        a_sum = 0
        for i in range(1, n):
            j = n - 1 - i
            if j > 0:
                a_sum += a[i] * m[j]
        a[n] = m[n] + a_sum
        
        # Calculate c(n)
        c_sum = 0
        for i in range(1, n):
            j = n - 1 - i
            if j > 0:
                c_sum += c[i] * a[j]
        c[n] = a[n] + c_sum

    print(f"Calculation for the number of expressions of length 5, c(5):")
    
    # Building the equation string for c(5)
    c1_a_val = c[1] * a[max_len - 1 - 1]
    c2_a_val = c[2] * a[max_len - 1 - 2]
    c3_a_val = c[3] * a[max_len - 1 - 3]
    # In general this would be a loop, but for n=5 we can be explicit
    # The sum part is for i from 1 to 3. i=4 gives j=0, which is excluded.
    # c(5) = a(5) + c(1)a(3) + c(2)a(2) + c(3)a(1)
    
    # Python code requires explicit access for values used in print statements.
    c1, c2, c3 = c[1], c[2], c[3]
    a1, a2, a3 = a[1], a[2], a[3]
    a5 = a[5]

    c1a3 = c1 * a3
    c2a2 = c2 * a2
    c3a1 = c3 * a1
    
    print(f"c(5) = a(5) + c(1)*a(3) + c(2)*a(2) + c(3)*a(1)")
    print(f"c(5) = {a5} + {c1}*{a3} + {c2}*{a2} + {c3}*{a1}")
    print(f"c(5) = {a5} + {c1a3} + {c2a2} + {c3a1}")
    
    final_result = a5 + c1a3 + c2a2 + c3a1
    print(f"c(5) = {final_result}")
    
    # To conform with the final output format.
    print(f"\nThe total number of true boolean expressions is {final_result}.")


solve_boolean_expressions()
<<<90>>>