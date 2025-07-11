import sympy as sp

def solve_separatrix():
    """
    This function analyzes the given system of differential equations to find the separatrix.
    The analysis is explained in the thought process and comments.
    """
    
    # Define variables
    u, d = sp.symbols('u d')
    
    # System of ODEs
    # u_prime = u**2 * (u - 1)
    # d_prime = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    
    # Step 1 & 2: Find and classify equilibrium points.
    # From u_prime = 0, we get u=0 or u=1.
    # If u=0, d_prime = 2*d**2 = 0 -> d=0. Equilibrium point: (0, 0).
    # If u=1, d_prime = 2*d**2 + 2*d = 0 -> d=0 or d=-1. Equilibrium points: (1, 0), (1, -1).
    # The Jacobian analysis shows that (1, -1) is a saddle point, which is key to finding a separatrix.

    # Step 3 & 4: Search for an algebraic invariant curve d = f(u).
    # We test a simple form d = c * u**p. Let's try d = c * u**2.
    # Substitute d = c * u**2 into the ODEs.
    # d/dt = (d(d)/du) * (du/dt) = (2*c*u) * (u**2 * (u - 1)) = 2*c*u**4 - 2*c*u**3
    # The RHS of the d' equation is:
    # 2*(c*u**2)**2 + (-3*u + 5*u**2)*(c*u**2) - u**3*(1 - u)
    # = 2*c**2*u**4 - 3*c*u**3 + 5*c*u**4 - u**3 + u**4
    # = (2*c**2 + 5*c + 1)*u**4 + (-3*c - 1)*u**3
    
    # Equating coefficients of the powers of u:
    # For u**4: 2*c == 2*c**2 + 5*c + 1  => 2*c**2 + 3*c + 1 = 0
    # For u**3: -2*c == -3*c - 1         => c = -1
    
    # We check if c = -1 satisfies the first equation:
    # 2*(-1)**2 + 3*(-1) + 1 = 2 - 3 + 1 = 0.
    # It does. So, an invariant curve is d = -1 * u**2.
    
    c = -1
    p = 2
    
    # Step 5: This curve d = -u**2 passes through the saddle point (1, -1).
    # Its slope at u=1 is -2, which matches the slope of the unstable manifold of the saddle point.
    # Therefore, d = -u**2 is the unstable manifold of the saddle point and is the principal separatrix.

    print("The separatrix of the system is an invariant curve that passes through the saddle point (1, -1).")
    print("This curve has the equation d = c * u**p.")
    print(f"By solving the system, we find the coefficients to be c = {c} and p = {p}.")
    print("\nThe final equation for the separatrix is:")
    
    # Final output as per instructions.
    # This format explicitly includes each number in the final equation.
    print(f"d = ({c}) * u**({p})")

solve_separatrix()
<<<d = -u**2>>>