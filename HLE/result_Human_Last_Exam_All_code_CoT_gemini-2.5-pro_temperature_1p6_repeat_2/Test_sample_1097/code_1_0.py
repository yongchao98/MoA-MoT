import sympy as sp

def define_set_m():
    """
    Derives and prints the definition of the set M required to prove the
    existence and uniqueness of the solution to the BVP u'' - exp(u) = 0, 
    u(0) = u(1) = 0, using the Banach Fixed-Point Theorem.
    """

    print("To solve the BVP using the Banach Fixed-Point Theorem, we first convert it into a fixed-point problem u = T(u).")
    print("The BVP u''(x) = exp(u(x)) with u(0)=u(1)=0 is equivalent to the integral equation:")
    print("  u(x) = -Integral[0 to 1] of G(x, s) * exp(u(s)) ds")
    print("where G(x, s) is the Green's function for the operator -d^2/dx^2 with u(0)=u(1)=0.")
    print("This defines the operator T(u)(x).")
    print("\nWe need to find a complete metric space M such that T: M -> M and T is a contraction on M.")
    print("Any solution u(x) must be non-positive, so we search for a solution in a set of the form:")
    print("  M = {u in C_0[0, 1] | -R <= u(x) <= 0}, for some R > 0.")
    print("\nStep 1: Ensure T maps M to M (T: M -> M).")
    print("For any u in M, we have u(x) <= 0, which means exp(u(x)) <= 1.")
    print("Therefore, ||T(u)||_inf <= max_x Integral[0 to 1] of G(x, s) ds.")
    
    # Define symbolic variables
    x = sp.symbols('x')

    # The integral of the Green's function, h(x) = Integral G(x, s) ds,
    # is the solution to -h''(x) = 1 with h(0)=h(1)=0.
    # This gives h(x) = (x - x^2) / 2.
    h_x = (x - x**2) / 2
    
    # Find the maximum value of h(x) on the interval [0, 1]
    h_prime = sp.diff(h_x, x)
    # The critical point is where the derivative is zero
    critical_points = sp.solve(h_prime, x)
    crit_point_val = critical_points[0]
    # The maximum value is h(x) evaluated at the critical point
    max_h = h_x.subs(x, crit_point_val)
    R = max_h

    print(f"\nThe integral of the Green's function is h(x) = {h_x}.")
    print(f"Its maximum value on [0, 1] is {max_h}, which occurs at x = {crit_point_val}.")
    print(f"Thus, ||T(u)||_inf <= {R}. To ensure T maps M to M, we must choose R >= {R}.")
    print(f"We will choose R = {R}.")

    print("\nStep 2: Show T is a contraction on this M.")
    print("The Lipschitz constant k for T is bounded by max(exp(c)) * max_x(Integral(G(x,s)ds)).")
    print("For u in M, exp(c) <= exp(0) = 1. So the contraction constant k is bounded by the same maximum value.")

    k = max_h
    print(f"The contraction constant k <= {k}.")
    print(f"Since k = {k} < 1, T is a contraction mapping on M.")

    print("\nConclusion: The set M is defined as follows:")
    
    # Final definition of the set M
    final_set_definition = f"M = {{u in C[0, 1] | u(0) = u(1) = 0 and -{R} <= u(x) <= 0 for x in [0, 1]}}"
    print(final_set_definition)


if __name__ == '__main__':
    define_set_m()