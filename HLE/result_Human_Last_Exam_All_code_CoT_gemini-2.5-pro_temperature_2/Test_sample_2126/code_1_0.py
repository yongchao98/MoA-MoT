import math

def solve_spacetime_property():
    """
    This function outlines the final calculation for the spacetime property.
    The derivation leads to the expression Q(t) = exp(-12t) + exp(-24t).
    The final value is found by evaluating this at t=0.
    """
    
    # At t=0, the exponential terms become exp(0), which is 1.
    term1 = math.exp(0)
    term2 = math.exp(0)
    
    # The final equation is the sum of these two terms.
    result = term1 + term2
    
    print("The final calculation for the quantity at t=0 is:")
    # We output each number in the final equation as requested.
    num1 = int(term1)
    num2 = int(term2)
    final_result = int(result)
    
    print(f"{num1} + {num2} = {final_result}")

solve_spacetime_property()