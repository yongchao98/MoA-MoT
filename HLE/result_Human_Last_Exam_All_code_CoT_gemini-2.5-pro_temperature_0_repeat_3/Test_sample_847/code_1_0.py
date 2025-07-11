import numpy as np
import math

def solve_and_calculate():
    """
    Solves the optimization problem to find the minimal cost per bit.
    """
    # The optimal probability 'p' for a "yes" answer is the real root of p^3 + p - 1 = 0.
    # The coefficients of the polynomial are [1, 0, 1, -1].
    coefficients = [1, 0, 1, -1]
    roots = np.roots(coefficients)
    
    # Find the single real root from the solutions.
    p = roots[np.isreal(roots)].real[0]
    
    # The minimal cost per bit for a Type 2 question is -log_p(2).
    # This can be calculated as -ln(2) / ln(p).
    cost_per_bit_type2 = -math.log(2) / math.log(p)
    
    # The cost per bit for a Type 1 question is 2.
    cost_per_bit_type1 = 2
    
    # The optimal strategy uses the question type with the minimum cost per bit.
    min_cost_per_bit = min(cost_per_bit_type1, cost_per_bit_type2)
    
    # The problem asks for the final answer, which is this minimal cost coefficient.
    # The total cost C(n) for large n is approximately min_cost_per_bit * n * log2(n).
    # We will print the components of this final equation.
    
    print("To find the minimal cost, we analyze the cost per bit of information.")
    print(f"For Type 1 questions (comparisons), the cost per bit is {cost_per_bit_type1}.")
    print("\nFor Type 2 questions (general), the optimal strategy is to ask questions")
    print("where the probability 'p' of a 'yes' answer is the root of the equation p^3 + p - 1 = 0.")
    print(f"The real root of this equation is p = {p:.4f}.")
    print(f"This gives a minimal cost per bit of -log_p(2) = {cost_per_bit_type2:.3f}.")
    
    print(f"\nSince {cost_per_bit_type2:.3f} < {cost_per_bit_type1}, the optimal strategy is to use Type 2 questions.")
    
    print("\nThe final equation for the minimal number of coins C(n) for large n is:")
    print(f"C(n) ≈ {min_cost_per_bit:.3f} * n * log2(n)")
    
    # The final answer is the coefficient, rounded to 3 decimal places.
    # This is handled by the print statement above.

solve_and_calculate()