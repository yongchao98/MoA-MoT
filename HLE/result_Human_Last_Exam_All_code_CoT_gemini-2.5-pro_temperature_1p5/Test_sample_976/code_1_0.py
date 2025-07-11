import sympy as sp

def solve_electromagnetic_problem():
    """
    This function presents the derived solution for the electric potential and field
    for a conductive sphere in a uniform electric field.

    The derivation involves solving Laplace's equation in two regions (inside and
    outside the sphere) and applying boundary conditions for steady-state current.
    The boundary conditions are:
    1. Continuity of potential at r=R: Phi_1 = Phi_2
    2. Continuity of normal current density at r=R: sigma_1 * d(Phi_1)/dr = sigma_2 * d(Phi_2)/dr

    Solving these leads to the following expressions.
    """

    # Using string formatting to present the final mathematical expressions.
    # These expressions are the result of the physical derivation described above.
    # The problem asks for the solution in the region outside the sphere (r > R).
    # We provide the solution for both regions as given in the answer choices.

    # Potential expressions
    phi_inside_str = r"-\frac{3 \sigma_2}{\sigma_1 + 2 \sigma_2} E_0 r \cos\theta"
    phi_outside_str = r"-E_0 \left( r - \frac{(\sigma_1 - \sigma_2) R^3}{(\sigma_1 + 2\sigma_2)r^2} \right) \cos\theta"

    # Electric field expressions
    E_inside_str = r"\frac{3\sigma_2 E_0}{\sigma_1 + 2\sigma_2} (\cos\theta \hat{r} - \sin\theta \hat{\theta})"
    E_outside_str = r"E_0 \left[ 1 + \frac{2(\sigma_1 - \sigma_2) R^3}{(\sigma_1 + 2\sigma_2) r^3} \right] \cos\theta \hat{r} - E_0 \left[ 1 - \frac{(\sigma_1 - \sigma_2) R^3}{(\sigma_1 + 2\sigma_2) r^3} \right] \sin\theta \hat{\theta}"

    print("The derived solution for the potential and electric field is as follows:")
    print("-" * 70)

    print("Electric Potential \u03A6(r, \u03B8):")
    print(f"For r < R (inside): \u03A6_1(r, \u03B8) = {phi_inside_str}")
    print(f"For r > R (outside): \u03A6_2(r, \u03B8) = {phi_outside_str}")

    print("\nElectric Field E(r, \u03B8):")
    print(f"For r < R (inside): E_1(r, \u03B8) = {E_inside_str}")
    print(f"For r > R (outside): E_2(r, \u03B8) =\n {E_outside_str}")
    print("-" * 70)
    print("These expressions match answer choice B.")

if __name__ == '__main__':
    solve_electromagnetic_problem()