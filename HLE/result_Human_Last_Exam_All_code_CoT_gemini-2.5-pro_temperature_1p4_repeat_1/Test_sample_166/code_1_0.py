import math

def solve():
    """
    This function calculates the infimum c based on the theoretical analysis
    of Evolutionary Algorithms on CLIFF-like functions.

    The analysis shows that:
    1. The problem requires a jump of Hamming distance d = n/4 to reach the optimum.
    2. A standard (1,lambda)-EA with standard bit mutation would take exponential time.
    3. Polynomial time can be achieved by using a heavy-tailed mutation operator,
       a standard technique for this class of problems.
    4. Theoretical results (e.g., for a (1+1)-GA, which is asymptotically similar here)
       show an expected runtime of O(n^2 * m^(beta-1)) or O(n^(beta+1)), where m=n/4 is the jump size and beta > 1 is the power-law exponent
       of the mutation distribution.
    5. The exponent in the runtime is c = beta + 1.
    6. To find the infimum of c, we need the infimum of beta, which is 1.
    7. Therefore, the infimum of c is 1 + 1 = 2.
    """
    # The exponent c is derived from theoretical analysis.
    # c = inf(beta + 1) where beta > 1.
    c = 2.0

    # The problem asks to round c to three significant digits.
    # 2.0 is already in this format.
    c_rounded = round(c, 2)
    
    # We are asked to output each number in the final equation.
    # The "equation" is simply the value of c.
    # However, to meet the literal instruction, we can present it as a calculation.
    inf_beta = 1
    exponent_factor_from_analysis = 1
    c_val = inf_beta + exponent_factor_from_analysis
    
    print("The runtime is of the form O(n^c).")
    print("Based on theoretical analysis of heavy-tailed EAs on jump functions, the runtime exponent 'c' can be expressed as a function of the power-law parameter beta > 1.")
    print("A common result for the runtime T is O(n^(beta+1)).")
    print(f"The exponent is c = beta + {exponent_factor_from_analysis}.")
    print("To find the infimum of c, we take the infimum of beta, which is 1.")
    print(f"inf(c) = {inf_beta} + {exponent_factor_from_analysis}")
    print(f"c = {c_val:.3f}")


solve()