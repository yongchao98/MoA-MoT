import sympy

def solve():
    """
    This script explains the process of finding the set M and uses sympy
    to calculate the contraction constant for the Banach Fixed-Point Theorem proof.
    """
    # Define symbols for our calculations
    x, s = sympy.symbols('x s')

    print("Step 1: Reformulate the BVP as a fixed-point problem u = T(u).")
    print("The BVP u''(x) = exp(u(x)) with u(0)=u(1)=0 is equivalent to the integral equation:")
    print("u(x) = T(u)(x) = Integral[0 to 1] G(x, s) * exp(u(s)) ds")
    print("where G(x,s) is the Green's function for the operator u'' with zero boundary conditions.")
    print("G(x,s) = (s-1)*x for x <= s")
    print("G(x,s) = s*(x-1) for s >= x")
    print("\nStep 2: Define the set M.")
    print("From u'' = exp(u) > 0, we know u(x) is convex. With u(0)=u(1)=0, it implies u(x) <= 0.")
    print("So, we choose M to be the set of non-positive continuous functions on [0,1].")
    print("M = {u in C[0,1] | u(x) <= 0 for all x in [0,1]}")
    print("M is a closed subset of the complete space C[0,1], so M is also a complete metric space.")

    print("\nStep 3: Verify the conditions for the Banach Fixed-Point Theorem.")
    print("a) T maps M to M: For u in M, u(s)<=0, so exp(u(s))>0. Since G(x,s)<=0, the integrand is non-positive, thus T(u)(x) <= 0. So T(u) is in M.")
    print("b) T is a contraction on M:")
    print("||T(u) - T(v)|| <= ||u - v|| * max_x(Integral[0 to 1] |G(x,s)| ds)")
    print("The Lipschitz constant for exp(z) on z<=0 is 1. We need to calculate the contraction constant k.")

    # |G(x,s)| = s(1-x) for s <= x
    # |G(x,s)| = x(1-s) for s >= x
    # Integrate |G(x,s)| w.r.t. s from 0 to 1, splitting at s=x
    integral_of_G = sympy.simplify(sympy.integrate(s * (1 - x), (s, 0, x)) + sympy.integrate(x * (1 - s), (s, x, 1)))

    print(f"\nThe integral of |G(x,s)| w.r.t. s from 0 to 1 is: {integral_of_G}")

    # The contraction constant k is the maximum of this function on [0,1]
    k = sympy.maximum(integral_of_G, x, sympy.Interval(0, 1))

    print(f"The maximum value of this integral, which is our contraction constant k, is: {k}")

    # Final conclusion based on the value of k
    print("\nSince k = 1/8 is less than 1, T is a contraction mapping on M.")
    print("The final contraction inequality is: ||T(u) - T(v)|| <= 1/8 * ||u - v||")
    print("\nBy the Banach Fixed-Point Theorem, there exists a unique global solution u in the set M.")
    print("\nThe correct set M is therefore the set of non-positive continuous functions on [0,1].")

solve()