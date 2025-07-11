import sys

def solve_tripled_fixed_point_conditions():
    """
    This script explains and prints the conditions for the existence and
    uniqueness of an FGH-tripled fixed point.
    """

    # Step 1: Define the FGH-tripled fixed point
    print("An FGH-tripled fixed point is a triplet (x, y, z), where x is in X, y is in Y, and z is in Z,")
    print("that simultaneously satisfies the following system of equations:")
    print("1. F(x, y, z) = x")
    print("2. G(y, x, y) = y")
    print("3. H(z, y, x) = z  (Note: Assuming the second G in the prompt is a distinct function H)")
    print("-" * 50)

    # Step 2: State the main assumptions from Fixed-Point Theory
    print("The existence and uniqueness of such a point can be guaranteed by the Banach Fixed-Point Theorem.")
    print("The primary conditions are:")
    print("\nA. The sets X, Y, and Z must be complete metric spaces with metrics d_X, d_Y, and d_Z, respectively.\n")
    print("B. The functions F, G, and H must be Lipschitz continuous. This means they are 'well-behaved' and")
    print("their change in output is proportionally limited by the change in their inputs. Formally, they")
    print("must satisfy the following inequalities for some non-negative constants (alpha, beta, gamma):\n")

    # Step 3: Print the Lipschitz continuity conditions for F, G, and H
    print("For F: X*Y*Z -> X")
    print("d_X(F(x1, y1, z1), F(x2, y2, z2)) <= (alpha1 * d_X(x1, x2)) + (alpha2 * d_Y(y1, y2)) + (alpha3 * d_Z(z1, z2))")
    print("Numbers in equation: 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 2, 3, 1, 2\n")

    print("For G: Y*X*Y -> Y")
    print("d_Y(G(y1, x1, y1), G(y2, x2, y2)) <= (beta1 * d_Y(y1, y2)) + (beta2 * d_X(x1, x2)) + (beta3 * d_Y(y1, y2))")
    print("Numbers in equation: 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 2, 3, 1, 2\n")

    print("For H: Z*Y*X -> Z")
    print("d_Z(H(z1, y1, x1), H(z2, y2, x2)) <= (gamma1 * d_Z(z1, z2)) + (gamma2 * d_Y(y1, y2)) + (gamma3 * d_X(x1, x2))")
    print("Numbers in equation: 1, 1, 1, 2, 2, 2, 1, 1, 2, 2, 1, 2, 3, 1, 2")
    print("-" * 50)

    # Step 4: Print the final contraction conditions on the constants
    print("C. To guarantee that a combined operator T(x,y,z) = (F,G,H) is a contraction mapping,")
    print("the Lipschitz constants must collectively satisfy the following inequalities:\n")

    eq1_vars = ['alpha1', 'beta2', 'gamma3']
    eq2_vars = ['alpha2', 'beta1', 'beta3', 'gamma2']
    eq3_vars = ['alpha3', 'gamma1']

    # Final Condition 1
    print(f"Condition 1: {eq1_vars[0]} + {eq1_vars[1]} + {eq1_vars[2]} < 1")
    print("Numbers in equation: 1, 2, 3, 1")

    # Final Condition 2
    print(f"\nCondition 2: {eq2_vars[0]} + {eq2_vars[1]} + {eq2_vars[2]} + {eq2_vars[3]} < 1")
    print("Numbers in equation: 2, 1, 3, 2, 1")

    # Final Condition 3
    print(f"\nCondition 3: {eq3_vars[0]} + {eq3_vars[1]} < 1")
    print("Numbers in equation: 3, 1, 1")

if __name__ == '__main__':
    solve_tripled_fixed_point_conditions()