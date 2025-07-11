import math

def calculate_runtime_exponent():
    """
    Calculates the runtime exponent 'c' for the (1,lambda) EA on CLIFF_{3/4}.

    The analysis shows two phases for the algorithm:
    1. A 'ONEMAX'-like climb to the local optimum |x|_1 = 3n/4.
    2. A jump from the local optimum to the global optimum |x|_1 = n.

    The jump size 'd' is the number of bits that need to be flipped.
    d = n - 3n/4 = n/4.

    Standard mutation operators lead to superpolynomial runtime. However, advanced
    operators from "Fast EAs" can cross a cliff of width 'd' in O(n * d^1.5) evaluations.
    We use this result to find the runtime complexity.

    Complexity = O(n * (n/4)^1.5)
               = O(n * n^1.5)
               = O(n^(1 + 1.5))
               = O(n^2.5)

    The exponent 'c' is therefore 2.5.
    """

    # The runtime is dominated by the cliff-crossing phase.
    # The complexity is given by O(n * d^k) where d is the jump size and k is an exponent from theory.
    
    # Exponent of n in the O(n * ...) part of the complexity
    n_exponent_base = 1.0
    
    # Exponent of d in the O(... * d^k) part of the complexity, from Fast EA theory
    d_exponent = 1.5
    
    # The jump size d is a linear function of n, d = n/4.
    # So, (d)^k contributes k to the final exponent of n.
    n_exponent_from_d = d_exponent
    
    # The final exponent 'c' is the sum of the exponents of n.
    c = n_exponent_base + n_exponent_from_d
    
    print("Step-by-step calculation of the runtime exponent 'c':")
    print(f"The required jump size 'd' is n - 3n/4 = n/4.")
    print("Using results from Fast Evolutionary Algorithms, the runtime to cross the cliff is O(n * d^k).")
    print(f"The theoretical value for the exponent k is: {d_exponent}")
    print(f"Substituting d = n/4, the runtime becomes O(n * (n/4)^{d_exponent}) = O(n^1 * n^{d_exponent}) = O(n^(1 + {d_exponent})).")
    print(f"The exponent 'c' is the sum of the exponents of n.")
    print(f"c = {n_exponent_base} + {n_exponent_from_d} = {c}")
    
    # Rounding to three significant digits as requested
    c_rounded = round(c, 2) 
    
    print(f"\nThe infimum c is {c_rounded:.3g}.")


calculate_runtime_exponent()

# The final answer is the value of c
final_c = 2.5
# <<<2.50>>>