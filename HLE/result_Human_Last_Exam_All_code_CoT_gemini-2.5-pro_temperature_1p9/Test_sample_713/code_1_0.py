import math

def solve_simplex_polytope_volume(d: int):
    """
    Calculates the range of possible volumes for the polytope P constructed from a d-simplex T.

    The problem states that a d-simplex T has n=d(d-1)/2 edges. This is a typo.
    A d-simplex has d+1 vertices and (d+1)d/2 edges. We proceed with the correct
    number of edges.

    Args:
        d: The dimension of the simplex.
    """
    if not isinstance(d, int) or d < 1:
        print("Dimension d must be a positive integer.")
        return

    # Calculate the lower bound coefficient m_d = d!
    # This is achieved for an orthoscheme.
    lower_bound_coeff = math.factorial(d)

    # Calculate the upper bound coefficient M_d = d! * C(d-1, floor((d-1)/2))
    # This is achieved for a regular simplex.
    # C(n, k) is the binomial coefficient "n choose k".
    if d == 1:
        # For d=1, binom(0,0)=1
        upper_bound_coeff = lower_bound_coeff * 1
    else:
        n_for_binom = d - 1
        k_for_binom = (d - 1) // 2
        binom_val = math.comb(n_for_binom, k_for_binom)
        upper_bound_coeff = lower_bound_coeff * binom_val

    # "equation" from the prompt is interpreted as the final expression for the set
    print(f"The dimension of the simplex is d = {d}.")
    print(f"Let V be the volume of the simplex T.")
    
    if lower_bound_coeff == upper_bound_coeff:
        print(f"The volume of the polytope P is uniquely determined.")
        print(f"Final Equation: Vol(P) = {lower_bound_coeff} * V")
    else:
        print(f"The set of possible values for the volume of P is the interval [m_d * V, M_d * V].")
        print(f"The lower bound coefficient is m_d = d! = {d}! = {lower_bound_coeff}.")
        n = d-1
        k = (d-1)//2
        print(f"The upper bound coefficient is M_d = d! * C(d-1, floor((d-1)/2)) = {d}! * C({n}, {k}) = {upper_bound_coeff}.")
        print(f"Final Equation: Vol(P) in [{lower_bound_coeff}*V, {upper_bound_coeff}*V]")

# Example usage with a dimension d. You can change this value.
# The original problem does not specify a value for d. Let's use d=4 as an example.
d_example = 4
solve_simplex_polytope_volume(d_example)
