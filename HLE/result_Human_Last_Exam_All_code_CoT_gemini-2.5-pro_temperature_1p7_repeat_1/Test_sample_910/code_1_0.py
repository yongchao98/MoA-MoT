import sympy as sp

def display_solution():
    """
    This function displays the chosen solution for the liquid rise height (xi)
    and the required voltage (V0).
    """
    # Define the symbols used in the equations
    xi, s, epsilon_0, V_0, rho, g, gamma = sp.symbols('xi s varepsilon_0 V_0 rho g gamma')

    # Expression for the height xi from Option C
    # xi_expr = s * (epsilon_0 * V_0**2 / (2 * rho * g * s**3) - gamma / (rho * g * s))
    xi_expr_latex = r"   \xi = s \left( \frac{\varepsilon_0 V_0^2}{2 \rho g s^3} - \frac{\gamma}{\rho g s} \right)"

    # Expression for the voltage V0 from Option C
    # V0_expr = sp.sqrt(4 * rho * g * s**3 / epsilon_0) * sp.sqrt(1 + 2 * gamma * s / (rho * g))
    V0_expr_latex = r"   V_0 = \sqrt{\frac{4 \rho g s^3}{\varepsilon_0}} \left( 1 + \frac{2 \gamma s}{\rho g} \right)^{1/2}"

    stability_statement = (
        "The interface becomes unstable if the surface tension cannot counteract "
        "the electrostatic forces, leading to oscillatory behavior."
    )

    print("The expression for the height xi is:")
    print(xi_expr_latex)
    # Printing numbers in the equation: 2, 2, 3
    print("Numbers in the xi equation: 2, 2, 3")
    print("\nThe voltage V0 when the liquid rise is xi = s/2 is:")
    print(V0_expr_latex)
    # Printing numbers in the equation: 4, 3, 2, 1, 2
    print("Numbers in the V0 equation: 4, 3, 2, 1, 2")
    print("\nStability Discussion:")
    print(stability_statement)


if __name__ == "__main__":
    display_solution()