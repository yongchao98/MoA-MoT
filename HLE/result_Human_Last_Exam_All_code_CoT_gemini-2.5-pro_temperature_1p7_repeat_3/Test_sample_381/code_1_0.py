import math

def solve():
    """
    This function calculates the symbolic upper-bound for ||B * Q_{0, M}||_inf.
    
    The derivation proceeds as follows:
    1. Start with the norm inequality ||A||_inf <= sqrt(N) * ||A||_2.
       Therefore, ||B * Q||_inf <= sqrt(N) * ||B * Q||_2.
    2. The L2-norm ||B * Q||_2 can be shown to be bounded by a constant that is independent of N.
       This constant depends on the properties of the graph and the attention mechanism.
       It is proportional to the sum of the perturbations, sum(delta_t).
    3. The condition beta > 0 implies that sum(delta_t) is finite and proportional to (1 - beta).
    4. Combining these, the L2-norm ||B * Q||_2 is bounded by C' * (1 - beta) for some constant C'.
    5. The final infinity-norm bound is thus C * (1 - beta) * sqrt(N), where C is a constant.
    """
    
    # The problem asks for the upper-bound expression.
    # Since the constant 'C', 'beta', and 'N' are symbols, we will print the symbolic formula.
    
    # Representing the components of the final formula symbolically.
    constant_part = "C"
    beta_part = "1 - beta"
    n_part = "sqrt(N)"
    
    # Constructing the final expression string.
    # The question is: what is the upper-bound for ||B * Q_{0, M}||_inf expressed as a factor of sqrt(N)?
    # The bound is of the form Constant * (1 - beta) * sqrt(N).
    # Since beta is the limit of beta_k, we use beta in the final expression for the limit M -> infinity.
    final_bound_expression = f"{constant_part} * ({beta_part}) * {n_part}"
    
    print(final_bound_expression)

solve()
<<<C * (1 - beta) * sqrt(N)>>>