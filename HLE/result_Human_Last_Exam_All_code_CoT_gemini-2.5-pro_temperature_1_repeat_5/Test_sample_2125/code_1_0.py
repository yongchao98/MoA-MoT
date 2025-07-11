import numpy as np
from scipy.optimize import root_scalar

def find_alpha_zero():
    """
    Solves for the value of alpha_0 where F(alpha_0) = 0.
    This simplifies to finding the positive root of the polynomial equation
    2*alpha^4 - 2*alpha^3 - 1 = 0.
    """
    
    # Define the polynomial equation f(alpha) = 0
    # 2*alpha^4 - 2*alpha^3 - 1 = 0
    c4, c3, c2, c1, c0 = 2, -2, 0, 0, -1
    
    # Define the function for the root finder
    f = lambda alpha: c4 * alpha**4 + c3 * alpha**3 + c2 * alpha**2 + c1 * alpha + c0
    
    # Print the equation whose root we are finding
    print("The final equation for alpha is:")
    print(f"{c4}*alpha^4 + {c3}*alpha^3 + {c2}*alpha^2 + {c1}*alpha + {c0} = 0")
    
    # We can analyze the function to find a bracket for the root.
    # f(1) = 2 - 2 - 1 = -1
    # f(2) = 2*16 - 2*8 - 1 = 32 - 16 - 1 = 15
    # A root exists between 1 and 2.
    bracket = [1.0, 2.0]
    
    # Use a numerical solver to find the root
    # root_scalar is a robust way to find roots of a single-variable function.
    # The 'brentq' method is efficient and guaranteed to converge if the signs
    # of the function at the bracket endpoints are different.
    sol = root_scalar(f, bracket=bracket, method='brentq')
    
    alpha_0 = sol.root
    
    print(f"\nThe largest value alpha_0 such that F(alpha_0) = 0 is the positive root of the equation.")
    print(f"The calculated value is: {alpha_0}")
    
    return alpha_0

if __name__ == '__main__':
    alpha_0 = find_alpha_zero()
    # The final answer should be returned in the specified format
    # print(f"\n<<<{alpha_0}>>>")
