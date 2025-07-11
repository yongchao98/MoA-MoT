import math

def print_energy_formulas():
    """
    Prints the symbolic formulas for time-averaged stored energy per unit area
    for the evanescent wave, based on the selected answer choice D.
    """

    # Define common symbolic parts of the formulas for clarity
    sin_theta_sq = "sin(theta)**2"
    n_sq = "n**2"
    common_sqrt = f"sqrt({n_sq}*{sin_theta_sq} - 1)"
    
    # Define the common denominator string from the problem's options
    denominator_str = (
        f"2 * (omega/c) * ({n_sq} - 1) * "
        f"((({n_sq} + 1)*{sin_theta_sq}) - 1) * {common_sqrt}"
    )

    # Numerator for Electric Field Energy from Option D
    numerator_E_field_str = f"{n_sq} * (2*{n_sq}*{sin_theta_sq} - 1)"
    
    # Numerator for Magnetic Field Energy from Option D
    numerator_H_field_str = f"{n_sq} * ({n_sq}*{sin_theta_sq} - 1)"

    # Final expressions as strings
    energy_E_str = f"({numerator_E_field_str}) / ({denominator_str}) * epsilon_0 * |E_x0_i|^2"
    energy_H_str = f"({numerator_H_field_str}) / ({denominator_str}) * epsilon_0 * |E_x0_i|^2"

    print("Based on option D, the formulas are:")
    print("\nEnergy in E field =")
    print(energy_E_str)
    
    print("\nEnergy in H field =")
    print(energy_H_str)

# Execute the function to print the formulas
print_energy_formulas()
