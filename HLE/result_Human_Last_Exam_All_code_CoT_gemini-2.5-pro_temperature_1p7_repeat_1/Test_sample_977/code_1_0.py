import numpy as np

def calculate_electric_potential(sigma_0, k, a, b, epsilon_1, epsilon_2, x, y):
    """
    Calculates the electric potential Phi(x, y) for the given configuration.

    The function first determines which region (y > 0 or y < 0) the point (x, y)
    is in and then applies the appropriate formula derived from electrostatic
    boundary value problem. It also prints the formula with the numerical
    values substituted in.
    """
    # Check for valid inputs
    if a <= 0 or b <= 0:
        print("Error: The distances 'a' and 'b' to the grounded plates must be positive.")
        return None

    # Common denominator term
    # Denominator = k * [e2*cosh(k*a)*sinh(k*b) + e1*sinh(k*a)*cosh(k*b)]
    common_denominator = k * (epsilon_2 * np.cosh(k * a) * np.sinh(k * b) +
                              epsilon_1 * np.sinh(k * a) * np.cosh(k * b))

    if common_denominator == 0:
        print("Error: Calculation resulted in division by zero. Check parameters.")
        return None

    # Region 2: 0 <= y <= a
    if 0 <= y <= a:
        # Numerator = -sigma_0 * sinh(k*b) * sinh(k*(y-a)) * sin(k*x)
        numerator = -sigma_0 * np.sinh(k * b) * np.sinh(k * (y - a)) * np.sin(k * x)
        potential = numerator / common_denominator
        
        print("The point (x,y) is in the region 0 <= y <= a.")
        print("The potential Phi(x, y) is calculated using the formula:")
        print("Phi = [-sigma_0 * sinh(k*b) * sinh(k*(y-a)) * sin(k*x)] / [k * (e2*cosh(k*a)*sinh(k*b) + e1*sinh(k*a)*cosh(k*b))]")
        print("\nSubstituting the given values:")
        # To avoid clutter, intermediate products are calculated first.
        # Numerator terms
        num_term1 = -sigma_0
        num_term2 = np.sinh(k*b)
        num_term3 = np.sinh(k*(y-a))
        num_term4 = np.sin(k*x)
        # Denominator terms
        den_term1 = k
        den_term2_1 = epsilon_2 * np.cosh(k*a) * np.sinh(k*b)
        den_term2_2 = epsilon_1 * np.sinh(k*a) * np.cosh(k*b)
        
        print(f"Phi = [{num_term1:.4f} * {num_term2:.4f} * {num_term3:.4f} * {num_term4:.4f}] / "
              f"[{den_term1:.4f} * ({epsilon_2:.4f}*{np.cosh(k*a):.4f}*{np.sinh(k*b):.4f} + "
              f"{epsilon_1:.4f}*{np.sinh(k*a):.4f}*{np.cosh(k*b):.4f})]")
        print(f"Phi = [{numerator:.4f}] / [{common_denominator:.4f}]")

    # Region 1: -b <= y < 0
    elif -b <= y < 0:
        # Numerator = sigma_0 * sinh(k*a) * sinh(k*(y+b)) * sin(k*x)
        numerator = sigma_0 * np.sinh(k * a) * np.sinh(k * (y + b)) * np.sin(k * x)
        potential = numerator / common_denominator

        print("The point (x,y) is in the region -b < y < 0.")
        print("The potential Phi(x, y) is calculated using the formula:")
        print("Phi = [sigma_0 * sinh(k*a) * sinh(k*(y+b)) * sin(k*x)] / [k * (e2*cosh(k*a)*sinh(k*b) + e1*sinh(k*a)*cosh(k*b))]")
        print("\nSubstituting the given values:")
        # Numerator terms
        num_term1 = sigma_0
        num_term2 = np.sinh(k*a)
        num_term3 = np.sinh(k*(y+b))
        num_term4 = np.sin(k*x)
        
        print(f"Phi = [{num_term1:.4f} * {num_term2:.4f} * {num_term3:.4f} * {num_term4:.4f}] / "
              f"[{k:.4f} * ({epsilon_2:.4f}*{np.cosh(k*a):.4f}*{np.sinh(k*b):.4f} + "
              f"{epsilon_1:.4f}*{np.sinh(k*a):.4f}*{np.cosh(k*b):.4f})]")
        print(f"Phi = [{numerator:.4f}] / [{common_denominator:.4f}]")

    else:
        print(f"Error: The point y={y} is outside the defined regions [-b, a]=[-{b}, {a}].")
        return None
    
    print("\n---")
    print(f"The calculated electric potential at (x={x}, y={y}) is: {potential:.6f} V")
    return potential

if __name__ == '__main__':
    # Example Parameters based on the problem statement
    # You can change these values to see how the potential changes
    sigma_0_val = 1.0       # Surface charge density amplitude (C/m^2)
    k_val = 1.0             # Wave number (rad/m)
    a_val = 1.0             # Distance to upper plate (m)
    b_val = 1.0             # Distance to lower plate (m)
    epsilon_1_val = 8.854e-12 * 1 # Permittivity of region 1 (F/m), e.g., vacuum
    epsilon_2_val = 8.854e-12 * 2 # Permittivity of region 2 (F/m), e.g., dielectric const 2

    # Point of interest (x, y)
    # The question asks for the region 0 <= y <= a
    x_val = np.pi / (2 * k_val)  # Choose x such that sin(kx) = 1 for max potential effect
    y_val = 0.5 * a_val          # A point within the region 0 < y < a

    print("Calculating potential for the region 0 <= y <= a")
    calculate_electric_potential(sigma_0_val, k_val, a_val, b_val,
                                 epsilon_1_val, epsilon_2_val, x_val, y_val)
    
    print("\n" + "="*50 + "\n")

    # Another example for the other region
    y_val_region1 = -0.5 * b_val
    print("Calculating potential for the region -b < y < 0")
    calculate_electric_potential(sigma_0_val, k_val, a_val, b_val,
                                 epsilon_1_val, epsilon_2_val, x_val, y_val_region1)