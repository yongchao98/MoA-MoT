import math

def solve_problem():
    """
    Solves the problem by deriving the general formula and then calculating the specific case.
    
    Part (a): Derivation of the maximum number of roots.
    Let W = W(f_1, ..., f_t). By choosing the exponents k_i and l_i to be
    sufficiently large and distinct integers (e.g., k_i >= t-1, l_i >= t-1),
    we ensure W(x) is a polynomial for which we can analyze the roots.
    
    1. The degree of the Wronskian polynomial W(x) is given by:
       deg(W) = sum(k_i + l_i) - t*(t-1)/2
    
    2. The multiplicity of the root at x=0 is:
       N_0 = sum(k_i) - t*(t-1)/2
       
    3. The multiplicity of the root at x=1 is:
       N_1 = sum(l_i) - t*(t-1)/2
       
    The number of roots in the open interval ]0, 1[ is the total number of roots (degree)
    minus the roots at the endpoints (N_0 + N_1).
    
    Max roots in ]0, 1[ = deg(W) - N_0 - N_1
                       = (sum(k_i + l_i) - t*(t-1)/2) - (sum(k_i) - t*(t-1)/2) - (sum(l_i) - t*(t-1)/2)
                       = sum(k_i) + sum(l_i) - t*(t-1)/2 - sum(k_i) + t*(t-1)/2 - sum(l_i) + t*(t-1)/2
                       = t*(t-1)/2
                       
    So for part (a), the expression is t*(t-1)/2.
    
    Part (b): Calculation for t = 5.
    We substitute t=5 into the formula.
    """
    
    # Part (a)
    a_expression = "t*(t-1)/2"

    # Part (b)
    t = 5
    t_minus_1 = t - 1
    numerator = t * t_minus_1
    b_result = numerator // 2

    # Print the final answer in the required format
    # "output each number in the final equation!"
    print(f"(a) {a_expression}; (b) {t} * ({t} - 1) / 2 = {b_result}")

solve_problem()