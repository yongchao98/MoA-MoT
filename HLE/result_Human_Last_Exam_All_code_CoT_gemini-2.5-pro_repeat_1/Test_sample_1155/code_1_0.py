import sympy as sp

def solve_normal_cone():
    """
    This function calculates the explicit representation of the normal cone for the given problem.
    """
    # Step 1: Define the problem variables and functions
    x1, x2, x3 = sp.symbols('x1 x2 x3')
    x = sp.Matrix([x1, x2, x3])
    x_star = sp.Matrix([2, 0, -1])

    g = sp.Matrix([
        (x1 - 1)**2 + x2**2 - 1,
        (x1 - 3)**2 + x2**2 - 1,
        x3 + 1,
        -x3 - 2
    ])
    n_ineq = len(g)

    print("Problem Setup:")
    print(f"Feasible set F is defined by g(x) <= 0 for x in R^3.")
    print(f"Point x^* = {x_star.T}")
    print("-" * 30)

    # Step 2: Identify active constraints at x^*
    print("Step 1: Identify active constraints at x^*")
    active_indices = []
    for i in range(n_ineq):
        val = g[i].subs([(x1, x_star[0]), (x2, x_star[1]), (x3, x_star[2])])
        status = '(active)' if val == 0 else '(inactive)'
        print(f"g_{i+1}(x^*) = (x_{1} - {1 if i==0 else 3})^2 + x_{2}^2 - 1 = {val} {status}" if i < 2 else f"g_{i+1}(x^*) = ... = {val} {status}")
        if val == 0:
            active_indices.append(i + 1)

    print(f"\nThe set of active constraint indices is I(x^*) = {active_indices}")
    print("-" * 30)

    # Step 3: Compute gradients of active constraints
    print("Step 2: Compute gradients of the active constraints at x^*")
    active_gradients = []
    for i in active_indices:
        grad_gi = g[i-1].jacobian(x)
        grad_gi_at_x_star = grad_gi.subs([(x1, x_star[0]), (x2, x_star[1]), (x3, x_star[2])])
        active_gradients.append(grad_gi_at_x_star.T) # Store as column vector
        print(f"Gradient of g_{i} at x^*, grad(g_{i}(x^*)):")
        sp.pprint(grad_gi_at_x_star.T)
        print()

    print("-" * 30)

    # Step 4: Define the normal cone and derive its explicit representation
    print("Step 3: Characterize the normal cone T_F^°(x^*)")
    print("The normal cone is the conic hull of the active gradients:")
    print("T_F^°(x^*) = { s in R^3 | s = mu_1 * grad(g_1(x*)) + mu_2 * grad(g_2(x*)) + mu_3 * grad(g_3(x*)), for mu_1, mu_2, mu_3 >= 0 }")

    mu1, mu2, mu3 = sp.symbols('mu1 mu2 mu3', nonnegative=True)
    s_expr = mu1 * active_gradients[0] + mu2 * active_gradients[1] + mu3 * active_gradients[2]

    print("\nAn element s = (s1, s2, s3) in the normal cone is given by:")
    # Showing the equation with the numbers from the gradients
    print(f"s1 = {active_gradients[0][0]}*mu_1 + {active_gradients[1][0]}*mu_2 + {active_gradients[2][0]}*mu_3 = {s_expr[0]}")
    print(f"s2 = {active_gradients[0][1]}*mu_1 + {active_gradients[1][1]}*mu_2 + {active_gradients[2][1]}*mu_3 = {s_expr[1]}")
    print(f"s3 = {active_gradients[0][2]}*mu_1 + {active_gradients[1][2]}*mu_2 + {active_gradients[2][2]}*mu_3 = {s_expr[2]}")

    print("\nFrom these equations, we deduce the conditions on s1, s2, s3:")
    print("1. For s2: The equation is s2 = 0. This is a strict equality.")
    print("2. For s3: The equation is s3 = mu_3. Since mu_3 >= 0, we have the condition s3 >= 0.")
    print("3. For s1: The equation is s1 = 2*mu_1 - 2*mu_2. Since mu_1 and mu_2 can be any non-negative real numbers, s1 can take any real value.")
    print("-" * 30)

    # Final explicit representation
    print("Final Answer: The explicit representation of the normal cone is:")
    print("T_F^°(x^*) = { s = (s1, s2, s3) in R^3 | s1 is any real number, s2 = 0, s3 >= 0 }")
    print("\nThe final equation for the components of a vector in the normal cone are:")
    print(f"s2 = 0")
    print(f"s3 >= 0")


if __name__ == '__main__':
    solve_normal_cone()