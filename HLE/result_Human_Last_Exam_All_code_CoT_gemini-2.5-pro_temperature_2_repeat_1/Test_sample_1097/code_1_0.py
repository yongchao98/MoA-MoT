import numpy as np

def describe_set_M():
    """
    Describes the set M required to prove the existence and uniqueness of solutions
    for the given boundary value problem using the Banach Fixed-Point Theorem.
    """

    print("To apply the Banach Fixed-Point Theorem to the boundary value problem:")
    print("u''(x) - exp(u(x)) = 0, with u(0) = u(1) = 0\n")

    print("Step 1: Reformulate the problem as a fixed-point equation u = T(u).")
    print("The BVP can be transformed into an integral equation using Green's function G(x, s):")
    print("u(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds")
    print("This defines the operator T as: (Tu)(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds\n")

    print("Step 2: Define the set M.")
    print("M must be a complete metric space where T is a contraction that maps M to itself.")
    print("We start with the space C_0[0,1] = {u in C[0,1] | u(0) = 0 and u(1) = 0}, equipped with the sup-norm ||u||_inf.")
    print("We define M as a closed ball of radius R within this space:")
    print("M = {u in C_0[0,1] | ||u||_inf <= R}")
    print("where R is a positive constant that we need to determine.\n")

    print("Step 3: Find conditions on R for T to be a contraction from M to M.")
    
    print("\n  Condition (a): T must map M to itself (i.e., for u in M, ||Tu||_inf <= R).")
    print("  Analysis shows this requires: exp(R) * sup_x | integral_0^1 G(x,s) ds | <= R.")
    print("  The integral part evaluates to 1/8. So, the final condition is:")
    print("  exp(R) / 8 <= R")

    print("\n  Condition (b): T must be a contraction on M (||Tu - Tv|| <= k * ||u - v|| with k < 1).")
    print("  Analysis shows the contraction constant k = exp(R) / 8.")
    print("  For T to be a contraction, we need k < 1, which leads to:")
    print("  exp(R) / 8 < 1\n")

    ln_8 = np.log(8)

    print("--- Final Definition of M ---")
    print("The set M is a closed ball of radius R in C_0[0,1]:")
    print("M = {u in C[0,1] | u(0) = u(1) = 0, ||u||_inf <= R}")
    print("\nwhere R is any constant that satisfies BOTH of the following inequalities:")
    print(f"1. exp(R) / 8 <= R")
    print(f"2. R < ln(8)  (which means R < {ln_8:.4f})\n")

    print("--- Example ---")
    print("A concrete example of a valid radius is R = 1.")
    r_example = 1.0
    cond1_val = np.exp(r_example) / 8
    print(f"For R={r_example}, let's check the conditions:")
    print(f"  1. Is exp({r_example}) / 8 <= {r_example}?  ->  {cond1_val:.4f} <= {r_example}. This is true.")
    print(f"  2. Is {r_example} < ln(8)?         ->  {r_example} < {ln_8:.4f}. This is true.")
    print("\nTherefore, a suitable choice for the set is M = {u in C[0,1] | u(0) = u(1) = 0, ||u||_inf <= 1}.")

if __name__ == '__main__':
    describe_set_M()