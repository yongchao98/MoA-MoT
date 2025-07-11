import math

def print_energy_formulas():
    """
    This function prints the symbolic formulas for the time-averaged stored energy
    per unit area in the electric and magnetic fields of an evanescent wave for p-polarization.
    """

    # Common denominator part for the expressions
    denominator_str = "2 * (w/c) * (n**2 - 1) * ((n**2 + 1) * sin(theta)**2 - 1) * sqrt(n**2 * sin(theta)**2 - 1)"

    # Numerator for the electric field energy
    electric_numerator_str = "n**2 * (2*n**2 * sin(theta)**2 - 1)"

    # Numerator for the magnetic field energy
    magnetic_numerator_str = "n**2 * (n**2 * sin(theta)**2 - 1)"

    # Final expressions
    electric_energy_str = f"Energy in E field = ({electric_numerator_str}) / ({denominator_str}) * epsilon_0 * |E_x0_i|^2"
    magnetic_energy_str = f"Energy in H field = ({magnetic_numerator_str}) / ({denominator_str}) * epsilon_0 * |E_x0_i|^2"

    print("The total time-averaged stored energy per unit area is given by:")
    print("-" * 60)
    print(electric_energy_str)
    print("")
    print(magnetic_energy_str)
    print("-" * 60)
    print("\nWhere:")
    print("  n: refractive index of the core")
    print("  theta: incident angle")
    print("  w: angular frequency of the light")
    print("  c: speed of light in vacuum")
    print("  epsilon_0: permittivity of free space")
    print("  |E_x0_i|^2: squared magnitude of the x-component of the incident electric field amplitude")

if __name__ == "__main__":
    print_energy_formulas()