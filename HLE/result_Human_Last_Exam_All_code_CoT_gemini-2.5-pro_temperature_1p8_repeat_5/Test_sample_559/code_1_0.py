import sympy

def find_separatrix():
    """
    This function outlines the steps to find the separatrix for the given system of ODEs
    and prints the final equation.
    """
    
    # Define the variables
    u, d = sympy.symbols('u d')
    
    # Define the system of ODEs
    du_dt = u**2 * (u - 1)
    dd_dt = 2*d**2 + (-3*u + 5*u**2)*d - u**3*(1-u)
    
    print("Step 1 & 2: Finding and Classifying Equilibrium Points")
    print("-----------------------------------------------------")
    
    # Find equilibrium points by setting derivatives to zero
    # From du/dt = u^2(u-1) = 0, we get u = 0 or u = 1.
    
    # Case u=0
    # dd_dt = 2*d^2 = 0  => d = 0
    eq1 = (0, 0)
    
    # Case u=1
    # dd_dt = 2*d^2 + 2*d = 2d(d+1) = 0 => d = 0 or d = -1
    eq2 = (1, 0)
    eq3 = (1, -1)
    
    print(f"The equilibrium points are: E1={eq1}, E2={eq2}, E3={eq3}")
    
    # Find the Jacobian matrix to classify the points
    F = du_dt
    G = dd_dt
    
    J = sympy.Matrix([[F.diff(u), F.diff(d)], [G.diff(u), G.diff(d)]])
    
    # Eigenvalues at E2 = (1, 0)
    J_E2 = J.subs({u: 1, d: 0})
    eigenvals_E2 = list(J_E2.eigenvals().keys())
    print(f"The Jacobian at E2={eq2} has eigenvalues {eigenvals_E2}, so it is an unstable node.")

    # Eigenvalues at E3 = (1, -1)
    J_E3 = J.subs({u: 1, d: -1})
    eigenvals_E3 = list(J_E3.eigenvals().keys())
    print(f"The Jacobian at E3={eq3} has eigenvalues {eigenvals_E3}, so it is a saddle point.")
    print("The separatrix is the invariant manifold passing through this saddle point.")
    print("\n")

    print("Step 3 & 4: Finding the Invariant Curve")
    print("-----------------------------------------")
    print("We hypothesize an invariant curve of the form d = a*u^2 + b*u + c.")
    print("By substituting this into the equation dd/du = (dd/dt) / (du/dt) and solving for the coefficients, we find two possible curves.")
    
    # As derived in the thought process, the two invariant curves are:
    curve1 = "-u**2"
    curve2 = "u*(1-u)"
    print(f"The two invariant polynomial curves are: d = {curve1} and d = {curve2}")
    
    print("\nStep 5: Identifying the Separatrix")
    print("------------------------------------")
    print("The separatrix must pass through the saddle point E3 = (1, -1).")
    
    # Check which curve passes through (1, -1)
    u_saddle, d_saddle = 1, -1
    
    # Check curve 1: d = -u^2
    d_check1 = -u_saddle**2
    print(f"For the curve d = -u^2, at u={u_saddle}, d = {d_check1}. This matches the saddle point E3=({u_saddle},{d_saddle}).")

    # Check curve 2: d = u(1-u)
    d_check2 = u_saddle*(1-u_saddle)
    print(f"For the curve d = u(1-u), at u={u_saddle}, d = {d_check2}. This does not match the saddle point.")
    
    print("\nTherefore, the separatrix is the curve d = -u^2.")
    print("\nFinal Answer")
    print("------------")
    # Represent d = -u^2 in the form d = a*u^2 + b*u + c and print the numbers.
    a = -1
    b = 0
    c = 0
    print("The equation for the separatrix is of the form: d = a*u**2 + b*u + c")
    print(f"The coefficients are: a = {a}, b = {b}, c = {c}.")
    print(f"So, the final equation is: d = ({a}) * u**2 + ({b}) * u + ({c})")
    
if __name__ == "__main__":
    find_separatrix()