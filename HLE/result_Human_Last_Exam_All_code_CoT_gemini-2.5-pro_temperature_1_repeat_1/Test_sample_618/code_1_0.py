import sympy
from sympy import symbols, Function, exp, Eq, dsolve, integrate, solve

def find_x_transformation():
    """
    This function derives and prints the finite transformation for the x-variable
    for the given PDE using the Lie symmetry method.
    """
    # Step 1: Define symbolic variables
    # t, x are the independent variables.
    # epsilon is the group parameter.
    # k1 is the constant from the PDE.
    # c1, c2, c3 are arbitrary constants from solving the determining equations.
    t, x, epsilon, k1, c1, c2, c3 = symbols('t x varepsilon k_1 c_1 c_2 c_3', real=True)

    # The new transformed variables are functions of the group parameter epsilon.
    t_star = Function('t_star')(epsilon)
    x_star = Function('x_star')(epsilon)

    print("Step 1: The infinitesimals for t and x are derived as:")
    # These are the results from solving the determining equations for the PDE.
    tau = c1
    xi = c2 * exp(k1 * t_star) + c3
    print(f"τ = {tau}")
    print(f"ξ = {c2 * exp(k1 * t) + c3}\n")


    print("Step 2: Set up the system of ordinary differential equations to find the finite transformations.")
    # The finite transformations are found by solving d(var*)/d(epsilon) = infinitesimal
    ode_t = Eq(t_star.diff(epsilon), tau)
    ode_x = Eq(x_star.diff(epsilon), xi)
    print("d(t*)/dε =", ode_t.rhs)
    print("d(x*)/dε =", ode_x.rhs, "\n")


    print("Step 3: Solve for the transformed variable x*.\n")

    # We need to distinguish between two cases for the constant c1.

    # Case 1: c1 is not equal to 0
    print("Case 1: c_1 ≠ 0")
    # Solve for t_star first
    t_star_sol_nonzero = dsolve(ode_t, t_star, ics={t_star.subs(epsilon, 0): t})
    t_star_expr = t_star_sol_nonzero.rhs

    # Substitute the solution for t_star into the ODE for x_star
    ode_x_subst_nonzero = ode_x.subs(t_star, t_star_expr)

    # Solve for x_star by integrating
    # The solution is x*(epsilon) = x*(0) + integral from 0 to epsilon
    integral_expr_nonzero = integrate(ode_x_subst_nonzero.rhs, (epsilon, 0, epsilon))
    x_star_sol_nonzero = x + integral_expr_nonzero

    print("The general transformation for x is:")
    # Using 'x_new' for the transformed variable for clarity in the output
    x_new = symbols('x\'')
    final_eq_nonzero = Eq(x_new, sympy.simplify(x_star_sol_nonzero))
    print(final_eq_nonzero)
    print(f"or, expanded:")
    print(f"{final_eq_nonzero.lhs} = {final_eq_nonzero.rhs.expand()}")
    print("-" * 20)

    # Case 2: c1 is equal to 0
    print("Case 2: c_1 = 0")
    # In this case, tau = 0, so t_star = t
    t_star_expr_zero = t
    
    # Substitute c1=0 into the ODE for x_star
    ode_x_subst_zero = ode_x.subs(c1, 0).subs(t_star, t_star_expr_zero)
    
    # Solve for x_star by integrating
    integral_expr_zero = integrate(ode_x_subst_zero.rhs, (epsilon, 0, epsilon))
    x_star_sol_zero = x + integral_expr_zero

    print("The general transformation for x is:")
    final_eq_zero = Eq(x_new, x_star_sol_zero)
    print(final_eq_zero)


if __name__ == '__main__':
    find_x_transformation()
    # The final answer is the derived formula for x'
    # We will format the non-zero case as the final answer string
    t, x, epsilon, k1, c1, c2, c3 = symbols('t x varepsilon k_1 c_1 c_2 c_3', real=True)
    x_new = symbols('x\'')
    rhs = x + c3*epsilon + c2*exp(k1*t)*(exp(k1*c1*epsilon) - 1)/(k1*c1)
    # The prompt requests the final answer in a special format.
    # Since the answer is a formula, we will provide the main case (c1 != 0)
    # The code prints both cases, which is more complete.
    final_answer_str = str(Eq(x_new, rhs)).replace("varepsilon", "ε")
    # print(f"\n<<<x' = x + c_3*ε + (c_2/(k_1*c_1))*exp(k_1*t)*(exp(k_1*c_1*ε) - 1)>>>")
    # Let's provide a slightly more readable form as the final final answer.
    # The code output itself is the primary answer. The final tag is a summary.
    print("\n<<<x' = x + c_3*ε + c_2*exp(k_1*t)*(exp(k_1*c_1*ε) - 1)/(k_1*c_1)>>>")
