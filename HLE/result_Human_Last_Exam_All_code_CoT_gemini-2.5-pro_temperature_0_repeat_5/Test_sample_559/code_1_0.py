def solve_separatrix():
    """
    This function describes and prints the equation of the separatrix for the given system.
    """
    print("The separatrix of the system is an invariant manifold passing through the saddle point (1, -1).")
    print("By solving the trajectory differential equation, we find the separatrix is a parabola.")
    print("The equation of the separatrix is of the form: d = a*u^2 + b*u + c")
    
    # Coefficients of the polynomial d = -u^2
    a = -1
    b = 0
    c = 0
    
    print("\nThe coefficients of the equation are:")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}")
    
    print("\nSo, the final equation for the separatrix is:")
    print(f"d = ({a}) * u^2 + ({b}) * u + ({c})")
    print("Which simplifies to: d = -u^2")

solve_separatrix()