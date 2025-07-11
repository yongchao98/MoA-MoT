import sympy

def find_separatrix():
    """
    This function outlines the analytical steps to find the separatrix of the given system
    and prints the final equation.
    """
    
    print("Step 1: Finding the equilibrium points.")
    print("The system of equations is:")
    print("d'(t) = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)")
    print("u'(t) = u**2*(u - 1)")
    print("Setting u'(t) = 0, we get u**2*(u - 1) = 0, which yields u = 0 or u = 1.")
    print("Case 1: u = 0. The first equation becomes d'(t) = 2*d**2 = 0, so d = 0. Equilibrium point: (0, 0).")
    print("Case 2: u = 1. The first equation becomes d'(t) = 2*d**2 + 2*d = 2*d*(d + 1) = 0, so d = 0 or d = -1. Equilibrium points: (0, 1) and (-1, 1).")
    print("The equilibrium points are (0, 0), (0, 1), and (-1, 1).\n")
    
    print("Step 2: Classifying the equilibrium points.")
    print("We analyze the stability by linearizing the system (using the Jacobian matrix).")
    print("This analysis shows that (0, 1) is an unstable node and (-1, 1) is a saddle point.")
    print("The point (0, 0) is a degenerate point.\n")

    print("Step 3: Identifying the separatrix.")
    print("The separatrix is a trajectory that separates regions of different behavior. The unstable manifold of the saddle point (-1, 1) acts as such a separatrix.\n")

    print("Step 4: Finding the equation for the separatrix.")
    print("We look for an invariant curve of the form d = f(u) that passes through the saddle point (-1, 1).")
    print("Let's test the hypothesis that the separatrix is a simple polynomial, d = -u**2.")
    print("To verify, we substitute d = -u**2 into the first differential equation:")
    print("The derivative d'(t) is, by the chain rule: d'(t) = -2*u * u'(t) = -2*u * (u**2*(u-1)) = -2*u**4 + 2*u**3.")
    print("The right-hand side of the first DE becomes: 2*(-u**2)**2 + (-3*u + 5*u**2)*(-u**2) - u**3*(1-u)")
    print("= 2*u**4 + 3*u**3 - 5*u**4 - u**3 + u**4 = -2*u**4 + 2*u**3.")
    print("Since the left-hand side equals the right-hand side, d = -u**2 is indeed a trajectory.\n")

    print("Step 5: Confirming it is the separatrix and providing the final equation.")
    print("We check if this curve passes through the saddle point (-1, 1).")
    print("At u = 1, d = -(1)**2 = -1. The curve passes through (-1, 1).")
    print("Therefore, the equation for the separatrix is d = -u**2.\n")
    
    print("Final Equation of the Separatrix:")
    # Define the components of the final equation
    # The equation is d = c * u^p
    c = -1
    p = 2
    
    print(f"The separatrix is given by the equation: d = {c}*u**{p}")
    print("\nThe numbers in this final equation are:")
    print(f"The coefficient of the u term (c): {c}")
    print(f"The power of the u term (p): {p}")

if __name__ == '__main__':
    find_separatrix()