import math

def calculate_theoretical_stress_at_notch():
    """
    Explains the concept of theoretical stress at a sharp notch tip.
    
    In linear elastic theory, the stress at the tip of a perfectly sharp notch
    (with a tip radius of zero) is singular, meaning it approaches infinity.
    This script explains this concept step by step.
    """
    
    # The nominal stress applied to the plate far from the notch.
    # We represent it as a string since it's a variable in the problem statement.
    sigma_nominal_str = "σ_y"
    
    # The maximum stress at the notch tip (point A).
    sigma_max_str = "σ_max"
    
    # The stress concentration factor.
    Kt_str = "K_t"
    
    # The tip radius of the notch. For a theoretically perfect sharp corner, it is zero.
    tip_radius_rho = 0
    
    print("This program explains the theoretical stress at the tip of a sharp wedge (point A).")
    print("-" * 70)
    
    print("Step 1: The relationship between maximum stress and nominal stress.")
    print(f"The maximum stress ({sigma_max_str}) at the tip of a notch is calculated using the nominal stress ({sigma_nominal_str}) and a stress concentration factor ({Kt_str}).")
    print(f"The formula is: {sigma_max_str} = {Kt_str} * {sigma_nominal_str}\n")
    
    print("Step 2: Determining the Stress Concentration Factor (K_t).")
    print(f"The factor {Kt_str} depends on the geometry of the notch, especially its tip radius (ρ).")
    print(f"The image shows a perfectly sharp corner, which means its theoretical tip radius is ρ = {tip_radius_rho}.\n")
    
    print("Step 3: Conclusion for a sharp notch.")
    print(f"According to linear elasticity theory, as the tip radius ρ approaches 0, the stress concentration factor {Kt_str} approaches infinity.")
    
    # Using the infinity symbol for demonstration.
    Kt_value_str = "∞ (infinity)"
    
    print(f"Therefore, for ρ = {tip_radius_rho}, we have {Kt_str} → {Kt_value_str}.\n")
    
    print("Final Calculation:")
    print("Based on the formula, we can write the final equation conceptually:")
    print(f"{sigma_max_str} = ({Kt_str}) * {sigma_nominal_str}")
    print(f"{sigma_max_str} = ({Kt_value_str}) * {sigma_nominal_str}")
    print(f"The result is that {sigma_max_str} is theoretically infinite.")
    
if __name__ == "__main__":
    calculate_theoretical_stress_at_notch()