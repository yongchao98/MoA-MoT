import sympy
from sympy import series, O

def solve_critical_exponent_order():
    """
    This script symbolically derives the order of the first correction
    to the critical exponent nu in phi^4 theory using the epsilon-expansion framework.
    """
    # Define the symbols.
    # u is the coupling constant.
    # epsilon is the parameter 4 - d, where d is the spatial dimension.
    # B and C are positive constants that arise from loop integrals. Their
    # exact values are not needed to determine the order of the expansion.
    u, epsilon = sympy.symbols('u epsilon')
    B, C = sympy.symbols('B C', positive=True)

    print("Step 1: Define the one-loop beta function β(u).")
    # The beta function describes the running of the coupling constant u.
    # To one-loop, β(u) = -ε*u + B*u^2.
    # The non-trivial fixed point u* is found by setting β(u*) = 0.
    beta_u = -epsilon * u + B * u**2
    print(f"β(u) = {beta_u}")
    u_star = sympy.solve(beta_u, u)[1] # The non-trivial solution u* != 0
    print(f"The non-trivial fixed point is at u* = {u_star}, which is of order O(ε).")
    print("-" * 20)

    print("Step 2: Define the one-loop anomalous dimension γ_φ²(u).")
    # This function is related to the renormalization of the φ² operator.
    # The one-loop diagram for this gives a term linear in u.
    # γ_φ²(u) = C*u + O(u^2)
    gamma_phi_squared = C * u
    print(f"γ_φ²(u) = {gamma_phi_squared}")
    print("This shows the leading term for the anomalous dimension is of order u^1.")
    print("-" * 20)

    print("Step 3: Use the RG relation for the critical exponent ν.")
    # The exponent ν is related to γ_φ² via the scaling dimension of the
    # 'thermal' operator (related to φ²). The relation is 1/ν = 2 - γ_φ²(u).
    nu_inverse_expr = 2 - gamma_phi_squared
    print(f"The inverse of ν is given by: 1/ν(u) = {nu_inverse_expr}")
    print("-" * 20)

    print("Step 4: Express ν(u) as a function of u and find its series expansion around u=0.")
    nu_expr = 1 / nu_inverse_expr
    print(f"This gives ν(u) = {nu_expr}")
    
    # Expand ν(u) in a Taylor series around u=0 to see the contributions order by order.
    nu_series = nu_expr.series(u, 0, 3)
    print("\nSeries expansion of ν(u) around u=0:")
    # The final equation with all its numbers:
    print(f"ν(u) = {nu_series}")
    print("-" * 20)
    
    print("Step 5: Analyze the series expansion to find the order of the first correction.")
    # The series is ν(u) = 1/2 + (C/4)*u + O(u^2).
    mean_field_value = nu_series.coeff(u, 0)
    first_correction_term = nu_series.coeff(u, 1) * u

    print(f"The mean-field value of ν (the u^0 term) is {mean_field_value}.")
    print(f"The first non-vanishing contribution beyond the mean-field value is '{first_correction_term}'.")
    print("This term is linear in u, meaning it is of order 1 in the coupling constant.")
    
    final_order = 1
    print(f"\nConclusion: The specific order is {final_order}.")

if __name__ == '__main__':
    solve_critical_exponent_order()
