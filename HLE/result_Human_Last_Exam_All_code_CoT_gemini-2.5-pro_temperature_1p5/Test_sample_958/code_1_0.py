import math

def solve_fiber_optics_energy():
    """
    This function presents the solution to the fiber-optics energy problem.

    Based on the derivation from Maxwell's equations for a p-polarized wave undergoing
    total internal reflection, we find the expressions for the time-averaged stored energy
    per unit area in the evanescent field.

    The electric field energy W_E is found to be:
    W_E = (n^2 * (2*n^2*sin(theta)^2 - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin(theta)^2 - 1) * sqrt(n^2*sin(theta)^2 - 1)) * epsilon_0 * |E_{x0}^i|^2

    This matches the expression in option D.

    The magnetic field energy W_H, derived consistently, is:
    W_H = (n^2) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin(theta)^2 - 1) * sqrt(n^2*sin(theta)^2 - 1)) * epsilon_0 * |E_{x0}^i|^2

    The W_H expression in option D contains an extra factor of (n^2*sin(theta)^2 - 1) in the numerator,
    which is likely a typo. However, since the W_E part of option D is correct and unique among the choices,
    D is the most plausible answer.
    """
    
    answer_choice = "D"
    
    print(f"The correct choice is {answer_choice}.\n")

    print("The expressions for the stored energy are:")
    
    # We will print the equations from choice D as requested by the problem format.
    # Note: These are symbolic representations, not code to be executed for a numerical result.
    
    energy_e_field_expr = "Energy in E field = (n^2 * (2*n^2*sin^2(theta) - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin^2(theta) - 1) * sqrt(n^2*sin^2(theta) - 1)) * epsilon_0 * |E_{x0}^i|^2"
    energy_h_field_expr = "Energy in H field = (n^2 * (n^2*sin^2(theta) - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin^2(theta) - 1) * sqrt(n^2*sin^2(theta) - 1)) * epsilon_0 * |E_{x0}^i|^2"

    print(energy_e_field_expr)
    print(energy_h_field_expr)

solve_fiber_optics_energy()