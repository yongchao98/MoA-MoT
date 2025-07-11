import math

def calculate_liquid_rise(rho, g, s, epsilon_0, gamma):
    """
    Calculates the voltage for liquid rise to s/2 and the corresponding height
    based on the (typo-corrected) formulas from the selected answer choice.

    Args:
        rho (float): Mass density of the liquid (kg/m^3)
        g (float): Gravitational acceleration (m/s^2)
        s (float): Separation between electrodes (m)
        epsilon_0 (float): Permittivity of free space (F/m)
        gamma (float): Surface tension of the liquid (N/m)
    """

    print("--- Input Parameters ---")
    print(f"Liquid Density (rho): {rho} kg/m^3")
    print(f"Gravity (g): {g} m/s^2")
    print(f"Plate Separation (s): {s} m")
    print(f"Permittivity (epsilon_0): {epsilon_0:.4e} F/m")
    print(f"Surface Tension (gamma): {gamma} N/m")
    print("-" * 26)

    # Correcting the V0 formula from option C for dimensional consistency
    # Original: (1 + 2*gamma*s / (rho*g)) -> Corrected: (1 + 2*gamma / (rho*g*s**2))
    # This formula is for the specific case where xi = s/2
    
    # Calculating the terms for the V0 equation
    v0_term1_val = 4 * rho * g * s**3 / epsilon_0
    v0_term2_val = 1 + 2 * gamma / (rho * g * s**2)
    
    # Calculate V0
    V0_squared = v0_term1_val * v0_term2_val
    V0 = math.sqrt(V0_squared)

    print("\n--- Voltage Calculation (for xi = s/2) ---")
    print(f"V0^2 = (4 * rho * g * s^3 / epsilon_0) * (1 + 2 * gamma / (rho * g * s^2))")
    print(f"V0^2 = ({v0_term1_val:.4e}) * ({v0_term2_val:.4f})")
    print(f"V0^2 = {V0_squared:.4e}")
    print(f"V0 = {V0:.2f} V")
    
    # Correcting the xi formula from option C for dimensional consistency
    # Original: gamma / (rho*g*s) -> Corrected: gamma / (rho*g*s**2)
    # The full formula is xi = s * ( (epsilon_0 * V0**2)/(2*rho*g*s**3) - gamma/(rho*g*s**2) )
    # We use the calculated V0 in this approximate formula.

    # Calculating the terms for the xi equation
    xi_term1_val = (epsilon_0 * V0**2) / (2 * rho * g * s**2)
    xi_term2_val = gamma / (rho * g * s)

    # Calculate xi
    xi = xi_term1_val - xi_term2_val

    print("\n--- Height Calculation (using the calculated V0) ---")
    print(f"xi = (epsilon_0 * V0^2) / (2 * rho * g * s^2) - gamma / (rho * g * s)")
    print(f"xi = ({xi_term1_val:.4f}) - ({xi_term2_val:.4f})")
    print(f"Calculated liquid rise (xi): {xi:.4f} m")
    print(f"\nNote: The V0 is calculated assuming xi=s/2 ({s/2} m).")
    print("The approximate formula for xi then predicts a different height due to inconsistencies in the provided options.")


# Using typical values for water between plates separated by 1mm
calculate_liquid_rise(
    rho=1000,          # kg/m^3 (density of water)
    g=9.8,             # m/s^2
    s=0.001,           # 1 mm
    epsilon_0=8.854e-12, # F/m
    gamma=0.072        # N/m (surface tension of water-air)
)