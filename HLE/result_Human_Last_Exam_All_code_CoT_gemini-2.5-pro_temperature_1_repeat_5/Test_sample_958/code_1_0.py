import math

def print_energy_formulas():
    """
    This function prints the symbolic formulas for the time-averaged stored energy
    per unit area for the electric and magnetic fields of an evanescent wave.
    The formulas correspond to the derived answer.
    """

    # Define the components of the formulas as strings for clear output.
    # These are not calculations, but representations of the symbolic formula.
    
    # Common factors and terms
    n_sq = "n^2"
    sin_theta_sq = "sin^2(theta)"
    E_term = "epsilon_0 * |E_x0^i|^2"
    
    # Numerator for Electric Field Energy
    numerator_E_str = f"{n_sq} * (2*{n_sq}*{sin_theta_sq} - 1)"
    
    # Numerator for Magnetic Field Energy
    numerator_H_str = f"{n_sq} * ({n_sq}*{sin_theta_sq} - 1)"
    
    # Denominator (common for both)
    denominator_str = f"2*(omega/c)*(n^2 - 1)*[({n_sq} + 1)*{sin_theta_sq} - 1]*sqrt({n_sq}*{sin_theta_sq} - 1)"

    # Print the final equations
    print("The total time-averaged stored energy per unit area is given by:")
    print("\nEnergy in E field = ")
    print(f"  ({numerator_E_str})")
    print(f"  --------------------------------------------------------------------------------------- * {E_term}")
    print(f"  ({denominator_str})")

    print("\nEnergy in H field = ")
    print(f"  ({numerator_H_str})")
    print(f"  --------------------------------------------------------------------------------------- * {E_term}")
    print(f"  ({denominator_str})")

# Execute the function to display the results
print_energy_formulas()