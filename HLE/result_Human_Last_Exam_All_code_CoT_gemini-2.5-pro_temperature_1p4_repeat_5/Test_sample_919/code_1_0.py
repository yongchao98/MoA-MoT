import numpy as np

def calculate_emi_shielding_force():
    """
    Calculates the force per unit area on a conducting plane due to a current sheet,
    as described in the EMI shielding problem.
    """
    # --- 1. Define Parameters and Constants ---
    # User-definable parameters
    K0 = 10.0      # Amplitude of the current sheet (A/m)
    a = 0.1        # Spatial variation constant (1/m)
    y = 5.0        # Position along the y-axis (m)
    d = 0.01       # Thickness of the air gap (m)
    mu_r = 80000.0 # Relative permeability of the magnetic material (e.g., Mu-metal)

    # Physical constant
    mu_0 = 4 * np.pi * 1e-7 # Permeability of free space (H/m)

    # --- 2. Calculate Intermediate Values ---
    mu = mu_r * mu_0  # Permeability of the magnetic material
    ad = a * d
    ay = a * y

    # --- 3. Calculate Numerator and Denominator of the Force Formula ---
    # Numerator: mu_0/2 * K0^2 * sin^2(ay)
    numerator_val = (mu_0 / 2) * (K0**2) * (np.sin(ay)**2)

    # Denominator: [cosh(ad) + (mu_0/mu) * sinh(ad)]^2
    cosh_ad = np.cosh(ad)
    sinh_ad = np.sinh(ad)
    mu_ratio = mu_0 / mu
    denominator_val = (cosh_ad + mu_ratio * sinh_ad)**2

    # --- 4. Calculate the Final Force ---
    force_magnitude = numerator_val / denominator_val

    # --- 5. Print the Results ---
    print("This script calculates the force per unit y-z area on the x=d interface.")
    print("-" * 70)
    
    print("Derived Formula (from option C):")
    print("f/A = (μ₀/2 * K₀² * sin²(ay)) / (cosh(ad) + (μ₀/μ) * sinh(ad))² î_x")
    print("-" * 70)

    print("Parameters Used:")
    print(f"  K₀ = {K0} A/m")
    print(f"  a = {a} 1/m")
    print(f"  y = {y} m")
    print(f"  d = {d} m")
    print(f"  μᵣ = {mu_r}")
    print(f"  μ₀ = {mu_0:.4e} H/m")
    print(f"  μ = μᵣ * μ₀ = {mu:.4e} H/m")
    print("-" * 70)

    print("Equation with numerical values:")
    
    # Building the string for the equation with numbers
    num_str = f"({mu_0:.4e} / 2) * ({K0})² * sin²({a} * {y})"
    den_str = f"(cosh({a} * {d}) + ({mu_0:.4e} / {mu:.4e}) * sinh({a} * {d}))²"
    
    print(f"f/A = [ {num_str} ] / [ {den_str} ] î_x")
    
    # Showing the evaluated parts of the equation
    num_calc_str = f"({numerator_val:.4e})"
    den_calc_str = f"({cosh_ad:.4f} + {mu_ratio:.4e} * {sinh_ad:.4f})² = ({denominator_val:.4f})"
    
    print(f"f/A = [ {num_calc_str} ] / [ {den_calc_str} ] î_x")
    print("-" * 70)
    
    print("Final Result:")
    print(f"The force per unit area is {force_magnitude:.4e} N/m² in the +x direction.")


if __name__ == "__main__":
    calculate_emi_shielding_force()