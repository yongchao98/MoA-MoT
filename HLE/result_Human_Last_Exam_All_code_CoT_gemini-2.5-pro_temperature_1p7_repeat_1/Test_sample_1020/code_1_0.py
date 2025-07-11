import math

def calculate_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_{n,d}.

    The formula for the complexity C is C = 2 + sum_{k=1}^{d-1} C(n, k),
    where C(n, k) is the binomial coefficient "n choose k".

    Args:
      n: The number of available indices.
      d: The number of variables in each monomial.
    """
    if not (2 <= d <= n):
        print("Error: The condition 2 <= d <= n must be met.")
        return

    # Although the problem states d is odd, the formula works for any d.
    # We can add a check if needed.
    # if d % 2 == 0:
    #     print("Warning: d is not odd, but proceeding with calculation.")
        
    complexity_sum_terms = []
    
    # Build the equation string part
    equation_str = "Complexity = 2"
    for k in range(1, d):
        equation_str += f" + C({n},{k})"
        
    print(equation_str)
    print("-" * len(equation_str))
    
    total_sum = 0
    for k in range(1, d):
        term = math.comb(n, k)
        complexity_sum_terms.append(str(term))
        print(f"C({n},{k}) = {term}")
        total_sum += term
        
    final_complexity = 2 + total_sum
    
    if complexity_sum_terms:
        print(f"Total complexity = 2 + {' + '.join(complexity_sum_terms)} = {final_complexity}")
    else: # This happens if d=2, loop range(1,2) gives k=1 but loop range(1,1) is empty. Oh wait d>=2, so d-1 >=1. range(1,d) is fine.
        print(f"Total complexity = {final_complexity}")
        
# Example from the problem statement: n, d are integers with 2 <= d <= n and d odd.
# Let's use n = 10 and d = 5 as an example.
n_example = 10
d_example = 5
calculate_complexity(n_example, d_example)