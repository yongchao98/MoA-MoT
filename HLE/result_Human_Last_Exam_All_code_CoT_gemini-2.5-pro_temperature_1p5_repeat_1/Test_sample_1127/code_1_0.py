import numpy as np

def solve_task():
    """
    This function provides the minimal polynomial for the connective constant
    of a width-2 triangular lattice strip.
    """
    
    # Coefficients of the minimal polynomial P(x) = x^3 - 2x^2 - 2x - 2
    # The list contains coefficients for x^3, x^2, x^1, and the constant term x^0.
    coeffs = [1, -2, -2, -2]
    
    equation_parts = []
    degree = len(coeffs) - 1
    
    print("The minimal polynomial for the connective constant (x) is:")
    
    for i, coeff in enumerate(coeffs):
        power = degree - i
        
        # Skip terms with a zero coefficient
        if coeff == 0:
            continue
            
        # Format the coefficient part, including parentheses as requested
        part = f"({coeff})"
        
        # Format the variable part
        if power > 1:
            part += f"*x^{power}"
        elif power == 1:
            part += f"*x"
            
        equation_parts.append(part)
        
    # Join the parts with ' + ' and form the equation
    equation_str = " + ".join(equation_parts) + " = 0"
    
    print(equation_str)

    # We can also use numpy to find the numerical value of the connective constant
    p = np.poly1d(coeffs)
    roots = p.roots
    
    # The connective constant is the largest (and only) positive real root
    connective_constant = None
    for root in roots:
        if np.isreal(root) and root > 0:
            connective_constant = np.real(root)
            break
            
    if connective_constant is not None:
        print(f"\nThe numerical value of the connective constant is approximately: {connective_constant}")

solve_task()