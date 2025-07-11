import numpy as np

def calculate_chi(k, h):
    """
    Calculates the stable amplitude factor chi based on the derived formula.
    The factor chi is a function of the dimensionless parameter kh.
    chi(kh) = kh / cosh(kh)
    
    Args:
        k (float): The wavenumber.
        h (float): The thickness of the plate.
        
    Returns:
        float: The value of the stable amplitude factor chi.
    """
    kh = k * h
    # Handle the case where kh is very large to avoid overflow in cosh
    if kh > 700:
        # For large kh, cosh(kh) is approximately 0.5 * exp(kh)
        # chi becomes 2 * kh * exp(-kh), which approaches 0
        return 0.0
    
    chi = kh / np.cosh(kh)
    return chi

def main():
    """
    Main function to calculate and print the result.
    Since k and h are not specified, we will use a representative
    case where the dimensionless product kh = 1.
    """
    # Representative case: k*h = 1.0
    # We can use k=1 and h=1 for this calculation.
    k_val = 1.0
    h_val = 1.0
    
    chi_value = calculate_chi(k_val, h_val)
    
    # The final equation relates the surface topography amplitude e_s
    # to the basal shear stress amplitude S_0.
    # The equation is e_s = chi * (S_0 / (delta_rho * g))
    
    print("For the representative case where k*h = 1:")
    print(f"The stable amplitude factor chi is: {chi_value:.4f}")
    print("\nThe final equation is:")
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # The only number we calculated is chi.
    print(f"e_s = {chi_value:.4f} * (S_0 / (delta_rho * g))")

if __name__ == "__main__":
    main()
