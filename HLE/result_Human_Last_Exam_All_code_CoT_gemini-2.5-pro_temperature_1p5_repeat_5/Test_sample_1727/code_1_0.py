import numpy as np

def calculate_chi(kh):
    """
    Calculates the stable amplitude factor chi for a given dimensionless wavenumber-thickness product kh.
    
    The calculation is based on the analytical solution to the Stokes flow problem with the specified boundary conditions.
    An ambiguity in the problem's boundary conditions at z=h is resolved by enforcing zero vertical velocity u_z(h)=0,
    a standard kinematic condition for channel flow.

    Args:
        kh (float): The dimensionless parameter k*h, where k is the wavenumber and h is the plate thickness.

    Returns:
        float: The value of the stable amplitude factor chi.
    """
    H = kh
    
    # Numerator of the expression for chi
    # Numerator = 4*H*cosh(H) - 3*sinh(H) - 3*H*exp(-H)
    term1_num = 4 * H * np.cosh(H)
    term2_num = -3 * np.sinh(H)
    term3_num = -3 * H * np.exp(-H)
    numerator = term1_num + term2_num + term3_num
    
    # Denominator of the expression for chi
    # Denominator = 2 * (sinh(H)^2 + H)
    denominator = 2 * (np.sinh(H)**2 + H)
    
    if denominator == 0:
        return np.inf # Avoid division by zero, though for H>0, this won't happen.
        
    chi = numerator / denominator
    return chi

def main():
    """
    Main function to calculate and print the amplitude factor chi.
    """
    # The amplitude factor chi is a function of the dimensionless parameter H = kh.
    # The derived analytical expression is:
    chi_formula = "(4*H*cosh(H) - 3*sinh(H) - 3*H*exp(-H)) / (2*(sinh(H)^2 + H))"
    
    print("The stable amplitude factor chi is given by the following equation, where H = kh:")
    print(f"chi(H) = {chi_formula}\n")

    # Example calculation for a specific value of kh, e.g., kh = 1.0
    kh_example = 1.0
    chi_value = calculate_chi(kh_example)

    print(f"For an example value of H = kh = {kh_example}:")
    
    H = kh_example
    # Print the full equation with numbers for the example case
    numerator_val = 4 * H * np.cosh(H) - 3 * np.sinh(H) - 3 * H * np.exp(-H)
    denominator_val = 2 * (np.sinh(H)**2 + H)

    # Output each number in the final equation
    print(f"chi({H}) = (4*{H:.1f}*cosh({H:.1f}) - 3*sinh({H:.1f}) - 3*{H:.1f}*exp(-{H:.1f})) / (2*(sinh({H:.1f})^2 + {H:.1f}))")
    print(f"chi({H}) = (4*{H:.1f}*{np.cosh(H):.4f} - 3*{np.sinh(H):.4f} - 3*{H:.1f}*{np.exp(-H):.4f}) / (2*({np.sinh(H):.4f}^2 + {H:.1f}))")
    print(f"chi({H}) = ({4*H*np.cosh(H):.4f} - {3*np.sinh(H):.4f} - {3*H*np.exp(-H):.4f}) / (2*({np.sinh(H)**2:.4f} + {H:.1f}))")
    print(f"chi({H}) = {numerator_val:.4f} / {denominator_val:.4f}")
    print(f"chi({H}) = {chi_value:.4f}")

if __name__ == '__main__':
    main()