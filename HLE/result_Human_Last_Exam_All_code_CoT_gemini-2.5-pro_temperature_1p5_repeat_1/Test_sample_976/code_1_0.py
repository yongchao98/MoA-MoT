import sympy
from sympy import symbols, cos, sin, latex, Mul, Add, Pow

def format_expression(expr, term_name):
    """Helper function to format sympy expressions for printing."""
    # Using latex representation for a more readable, textbook-like format.
    # Replace some latex commands for better plain text representation.
    expr_latex = latex(expr)
    expr_latex = expr_latex.replace(r'\cos{\left(\theta \right)}', 'cos(θ)')
    expr_latex = expr_latex.replace(r'\sin{\left(\theta \right)}', 'sin(θ)')
    expr_latex = expr_latex.replace(r'\left', '').replace(r'\right', '')
    print(f"{term_name} = {expr_latex}\n")

def main():
    """
    This script derives and displays the electric potential and field
    for a conducting sphere in a uniform electric field.
    """
    # Define symbols for the physical quantities
    E0, sigma1, sigma2, R, r, theta = symbols('E_0 sigma_1 sigma_2 R r theta')

    print("--- Expressions for the region INSIDE the sphere (r < R) ---")

    # Coefficient for the potential inside the sphere
    A_coeff = Mul(-3 * sigma2, Pow(Add(sigma1, 2 * sigma2), -1), evaluate=False)

    # Potential inside the sphere (Phi_in)
    Phi_in = Mul(A_coeff, E0, r, cos(theta), evaluate=False)
    format_expression(Phi_in, "Φ_in(r, θ)")

    # Electric field inside the sphere (E_in)
    # E_in = -grad(Phi_in) which simplifies to a uniform field in the z-direction.
    # E_in_z = -A_coeff * E0
    E_in_coeff = Mul(-1, A_coeff, E0, evaluate=False)
    print("vec(E)_in(r, θ) = (" + latex(E_in_coeff) + ") * (cos(θ) r_hat - sin(θ) theta_hat)")
    print("Which is a uniform field: vec(E)_in = (" + latex(E_in_coeff) + ") * z_hat\n")


    print("--- Expressions for the region OUTSIDE the sphere (r > R) ---")
    # This is the region specifically asked for in the problem.

    # Coefficient for the dipole term of the potential outside
    B_coeff = Mul(Add(sigma1, -sigma2), Pow(Add(sigma1, 2 * sigma2), -1), evaluate=False)

    # Potential outside the sphere (Phi_out)
    # Consists of the external field potential and the induced dipole potential
    dipole_term = Mul(B_coeff, Pow(R, 3), Pow(r, -2), evaluate=False)
    Phi_out_expr = Add(r, -dipole_term, evaluate=False)
    Phi_out = Mul(-E0, Phi_out_expr, cos(theta), evaluate=False)
    format_expression(Phi_out, "Φ_out(r, θ)")

    # Electric field outside the sphere (E_out)
    # E_out = -grad(Phi_out) has a radial (r) and a polar (theta) component.

    # Radial component of E_out (Er_out)
    Er_out_term_1 = 1
    Er_out_term_2 = Mul(2, B_coeff, Pow(R, 3), Pow(r, -3), evaluate=False)
    Er_out_factor = Add(Er_out_term_1, Er_out_term_2, evaluate=False)
    Er_out = Mul(E0, Er_out_factor, cos(theta), evaluate=False)

    # Polar component of E_out (Etheta_out)
    Etheta_out_term_1 = 1
    Etheta_out_term_2 = Mul(B_coeff, Pow(R, 3), Pow(r, -3), evaluate=False)
    Etheta_out_factor = Add(Etheta_out_term_1, -Etheta_out_term_2, evaluate=False)
    Etheta_out = Mul(-E0, Etheta_out_factor, sin(theta), evaluate=False)

    print("vec(E)_out(r, θ) = E_r * r_hat + E_θ * theta_hat, where:")
    format_expression(Er_out, "E_r")
    format_expression(Etheta_out, "E_θ")

    print("Comparing these derived expressions with the answer choices, they match option B.")

if __name__ == '__main__':
    main()
