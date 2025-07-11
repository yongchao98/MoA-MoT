import numpy as np

def analyze_flow_cytometry_setup(H_ext):
    """
    Simulates the effect of an external magnetic field on particles and the sensor.

    Args:
        H_ext (float): The external magnetic field from the permanent magnet (A/m).

    Returns:
        tuple: A tuple containing particle magnetization, particle saturation status,
               stray field at the sensor, and sensor saturation status.
    """
    # System Parameters
    # Susceptibility of paramagnetic particles
    particle_chi = 0.1
    # Saturation magnetization of a 1-um particle (A/m)
    particle_M_sat = 3000.0
    # Proportionality constant linking particle magnetization to sensor stray field (T*m/A)
    k_stray_field = 2.0e-7
    # Saturation field of the Spin Valve sensor (Tesla)
    sv_sensor_B_sat = 0.0005  # 0.5 mT

    # 1. Calculate Particle Magnetization
    # M is the lesser of the linear response and the saturation value.
    potential_M = particle_chi * H_ext
    actual_M = min(potential_M, particle_M_sat)

    # 2. Check for particle saturation
    is_particle_saturated = actual_M == particle_M_sat

    # 3. Calculate Stray Field at the Sensor
    # B_stray is proportional to the particle's actual magnetization.
    B_stray = k_stray_field * actual_M

    # 4. Check for sensor saturation
    is_sensor_saturated = B_stray > sv_sensor_B_sat

    return actual_M, is_particle_saturated, B_stray, is_sensor_saturated, particle_M_sat, potential_M, sv_sensor_B_sat

# --- Main Simulation ---

# Scenario 1: Properly positioned magnet
H_proper = 20000.0  # 20 kA/m
M_p, p_sat_p, B_p, s_sat_p, M_sat_p, pot_M_p, B_sat_p = analyze_flow_cytometry_setup(H_proper)

print("--- Scenario 1: Properly Positioned Magnet ---")
print(f"External Field H_ext = {H_proper:.0f} A/m")
print(f"Particle Magnetization M = min({pot_M_p:.1f}, {M_sat_p:.1f}) = {M_p:.1f} A/m")
print(f"Particle Saturated: {'YES' if p_sat_p else 'NO'}")
print(f"Resulting Stray Field B_stray = {B_p * 1000:.2f} mT")
print(f"Sensor Saturated (B_stray > {B_sat_p*1000:.2f} mT): {'YES' if s_sat_p else 'NO'}\n")


# Scenario 2: Improperly positioned (too strong) magnet
H_improper = 80000.0  # 80 kA/m
M_i, p_sat_i, B_i, s_sat_i, M_sat_i, pot_M_i, B_sat_i = analyze_flow_cytometry_setup(H_improper)

print("--- Scenario 2: Improperly Positioned Magnet ---")
print(f"External Field H_ext = {H_improper:.0f} A/m")
# Demonstrating Effect A: Particle Saturation
print(f"Particle Magnetization M = min({pot_M_i:.1f}, {M_sat_i:.1f}) = {M_i:.1f} A/m")
print(f"Particle Saturated: {'YES' if p_sat_i else 'NO'}")
# Demonstrating Effect D: Spin Valve Saturation
print(f"Resulting Stray Field B_stray = {B_i * 1000:.2f} mT")
print(f"Sensor Saturated (B_stray > {B_sat_i*1000:.2f} mT): {'YES' if s_sat_i else 'NO'}\n")

print("Conclusion:")
print("An improperly positioned (too strong) magnet first causes the paramagnetic particles to saturate (Effect A).")
print("These saturated particles then generate a stray field that is strong enough to saturate the highly sensitive Spin Valve sensor (Effect D).")
print("Sensor saturation is the most critical failure mode, as it renders the measurement data useless.")
print("Therefore, Spin Valve Saturation is the most significant negative effect.")
