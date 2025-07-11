import numpy as np

def calculate_chi(H):
    """
    Calculates the stable amplitude factor chi for a given dimensionless
    thickness H = kh.

    The formula is derived from solving the Stokes flow equations for a
    viscous layer with the specified boundary conditions.
    
    chi = (H * cosh(H)) / (cosh(H)^2 + 2*H*sinh(H)*cosh(H) - H^2)

    Args:
        H (float): The dimensionless thickness of the plate (H = kh).

    Returns:
        float: The value of the stable amplitude factor chi.
    """
    if H < 0:
        raise ValueError("H must be a non-negative number.")
    if H == 0:
        return 0.0

    cosh_H = np.cosh(H)
    sinh_H = np.sinh(H)
    
    numerator = H * cosh_H
    denominator = cosh_H**2 + 2 * H * sinh_H * cosh_H - H**2
    
    if denominator == 0:
        return np.inf
        
    chi = numerator / denominator
    return chi

def main():
    """
    Main function to get user input and print the result.
    An example calculation for H=1.0 is shown.
    """
    print("This script calculates the stable amplitude factor chi for a viscous plate model.")
    print("The factor chi relates surface deflection (e_s) to basal shear stress (S_0)")
    print("via the formula: e_s = chi * (S_0 / (delta_rho * g)).")
    
    try:
        # Example calculation for H=1.0
        H_example = 1.0
        print(f"\n--- Example Calculation for H = {H_example} ---")

        cosh_H_val = np.cosh(H_example)
        sinh_H_val = np.sinh(H_example)
        
        numerator_val = H_example * cosh_H_val
        denominator_val = cosh_H_val**2 + 2 * H_example * sinh_H_val * cosh_H_val - H_example**2
        
        chi_val = calculate_chi(H_example)

        print("\nThe formula for chi is: chi = (H * cosh(H)) / (cosh(H)^2 + 2*H*sinh(H)*cosh(H) - H^2)")
        print("\nPlugging in the numbers for H = 1.0:")
        
        # Outputting each number in the final equation as requested
        print(f"chi = ({H_example:.1f} * {cosh_H_val:.4f}) / (({cosh_H_val:.4f})^2 + 2*{H_example:.1f}*{sinh_H_val:.4f}*{cosh_H_val:.4f} - ({H_example:.1f})^2)")
        print(f"chi = ({numerator_val:.4f}) / ({denominator_val:.4f})")
        print(f"chi = {chi_val:.4f}")
        
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
