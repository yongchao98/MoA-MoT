import math

def calculate_weight_change(d, h, H, rho, t, g=9.8):
    """
    Calculates the estimated change in weight of a running hourglass based on a model
    where the impact force of the sand is the dominant effect.

    The formula used is: ΔW = (π * d^2 * h * ρ / (4 * t)) * sqrt(2 * g * (H - h))
    """
    
    # Convert units to SI (meters, kilograms, seconds)
    d_m = d / 100.0  # cm to m
    h_m = h / 200.0  # cm to m
    H_m = H / 100.0  # cm to m
    t_s = t * 60.0   # minutes to s

    # Calculate the mass flow rate (m_dot)
    # m_dot = (total mass) / time
    # total_mass = rho * volume = rho * (pi * (d/2)^2 * h)
    m_dot = (math.pi * d_m**2 * h_m * rho) / (4 * t_s)
    
    # Calculate the impact velocity (v_f)
    # v_f = sqrt(2 * g * L), where the characteristic fall distance L is taken as (H - h)
    fall_distance = H_m - h_m
    if fall_distance < 0:
        print("Error: The height of the sand pile (h) cannot be greater than the chamber height (H).")
        return

    impact_velocity = math.sqrt(2 * g * fall_distance)
    
    # Calculate the change in weight (Delta W)
    # Delta_W = m_dot * v_f
    delta_W = m_dot * impact_velocity
    
    # Print the step-by-step formula with the values
    print("The change in weight, ΔW, is estimated by the formula:")
    print("ΔW = (π * d² * h * ρ / (4 * t)) * sqrt(2 * g * (H - h))")
    print("\nSubstituting the given values (in SI units):")
    print(f"d = {d_m:.2f} m")
    print(f"h = {h_m:.2f} m")
    print(f"H = {H_m:.2f} m")
    print(f"ρ = {rho} kg/m³")
    print(f"t = {t_s} s")
    print(f"g = {g} m/s²")
    
    print("\nCalculation:")
    print(f"ΔW = (π * ({d_m:.2f})² * {h_m:.2f} * {rho} / (4 * {t_s})) * sqrt(2 * {g} * ({H_m:.2f} - {h_m:.2f}))")
    
    # Show intermediate calculations
    term1 = f"(3.14159 * {d_m**2:.4f} * {h_m:.2f} * {rho} / {4*t_s:.1f})"
    term2 = f"sqrt({2*g*fall_distance:.4f})"
    print(f"ΔW = ({m_dot:.6e}) * ({impact_velocity:.4f})")
    
    print(f"\nFinal Result:")
    print(f"The estimated change in weight ΔW is: {delta_W:.3e} N")
    print(f"This positive value indicates the hourglass is slightly heavier while running.")

# Given approximate values
d_cm = 1
h_cm = 2
H_cm = 4
rho_kg_m3 = 1500
t_min = 1

calculate_weight_change(d_cm, h_cm, H_cm, rho_kg_m3, t_min)
