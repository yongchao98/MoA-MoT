import numpy as np

def calculate_and_display_H_field(M0, Rp, R):
    """
    Calculates and displays the magnetic field H for the given setup.

    Args:
        M0 (float): Magnetization strength (in A/m).
        Rp (float): Radius of the magnetized sphere (in meters).
        R (float): Radius of the outer conductor cavity (in meters).
    """
    # First, display the final analytical equations corresponding to Option B.
    print("The correct expressions for the magnetic field H(r, θ) are:")
    print("=" * 70)
    
    # Region 1: Inside the shield (0 < r < Rp)
    # The equation is H = M0 * (2*Rp^3 + R^3)/(3*R^3) * (-cos(theta) i_r + sin(theta) i_theta)
    print(f"In the region 0 < r < {Rp}:")
    # Plugging the numbers into the equation string
    # Calculate the coefficient for H1
    h1_coeff_val = M0 * (2 * Rp**3 + R**3) / (3 * R**3)
    print(f"  H = {M0} * (2*{Rp}^3 + {R}^3) / (3*{R}^3) * [-cos(θ) î_r + sin(θ) î_θ]")
    print(f"  H = {h1_coeff_val:.4f} * [-cos(θ) î_r + sin(θ) î_θ] (A/m)")
    
    print("-" * 70)

    # Region 2: Between the shield and conductor (Rp < r < R)
    # The equations are:
    # H_r = - (2*M0/3) * [ (Rp/R)^3 - (Rp/r)^3 ] * cos(theta)
    # H_theta = (M0/3) * [ 2*(Rp/R)^3 + (Rp/r)^3 ] * sin(theta)
    print(f"In the region {Rp} < r < {R}:")
    # Plugging the numbers into the H_r equation string
    hr2_coeff_val = -2 * M0 / 3
    print(f"  H_r = - (2*{M0}/3) * [({Rp}/{R})^3 - ({Rp}/r)^3] * cos(θ)")
    print(f"  H_r = {hr2_coeff_val:.4f} * [{(Rp/R)**3:.4f} - ({Rp}/r)^3] * cos(θ) (A/m)")

    # Plugging the numbers into the H_theta equation string
    hth2_coeff_val = M0 / 3
    print(f"  H_θ = ({M0}/3) * [2*({Rp}/{R})^3 + ({Rp}/r)^3] * sin(θ)")
    print(f"  H_θ = {hth2_coeff_val:.4f} * [2*{(Rp/R)**3:.4f} + ({Rp}/r)^3] * sin(θ) (A/m)")
    print("=" * 70)

# --- Main execution ---
# Define the parameters for the problem
M0 = 1000000.0  # A/m (Example value for a strong permanent magnet)
Rp = 0.5         # meters
R = 1.0          # meters

# Call the function to display the solution with the given numbers
calculate_and_display_H_field(M0, Rp, R)