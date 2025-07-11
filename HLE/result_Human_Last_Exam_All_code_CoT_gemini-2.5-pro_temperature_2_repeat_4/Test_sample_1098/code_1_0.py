import math

# The apparatus described, with a one-way door, is a type of Brownian ratchet.
# It functions by rectifying the random thermal motion of gas particles.
# For particles to move from one compartment to the other, they must be in motion.

# The kinetic theory of gases tells us that the average kinetic energy of gas
# particles is directly proportional to the absolute temperature (in Kelvin).
# Particle speed is a direct consequence of this energy. A common measure is
# the root-mean-square (rms) speed, given by the formula:
# v_rms = sqrt(3 * k * T / m)
# where T is temperature, k is the Boltzmann constant, and m is the particle's mass.

# This script will calculate the rms speed for a nitrogen molecule at a non-zero
# temperature (room temperature) and at absolute zero to demonstrate that
# motion ceases at zero temperature, making temperature the crucial parameter.

# --- Constants ---
# Boltzmann constant in Joules/Kelvin
k = 1.38e-23
# Avogadro's number
N_A = 6.022e23
# Molar mass of Nitrogen (N2) in kg/mol
M_N2 = 28.0e-3
# Mass of a single nitrogen molecule in kg
m_N2 = M_N2 / N_A

# --- Scenario 1: Non-zero Temperature ---
T_room = 298.15  # Room temperature in Kelvin

# Calculation for T > 0 K
# v_rms^2 = (3 * k * T) / m
v_rms_sq_room = (3 * k * T_room) / m_N2
v_rms_room = math.sqrt(v_rms_sq_room)

print("To trap the gas, particles must move from one chamber to the other.")
print("This movement depends on their thermal energy, which is determined by temperature.\n")
print(f"Let's calculate the particle speed at a non-zero temperature (T = {T_room} K):")
print(f"v_rms = sqrt(3 * {k:.2e} J/K * {T_room} K / {m_N2:.2e} kg)")
print(f"Resulting speed: {v_rms_room:.2f} m/s")
print("At this temperature, particles move rapidly and can pass through the door.\n")

# --- Scenario 2: Absolute Zero Temperature ---
T_zero = 0  # Absolute zero in Kelvin

# Calculation for T = 0 K
# The formula directly shows that if T=0, the entire expression becomes 0.
v_rms_zero = 0

print(f"Now, let's calculate the particle speed at absolute zero (T = {T_zero} K):")
print(f"v_rms = sqrt(3 * {k:.2e} J/K * {T_zero} K / {m_N2:.2e} kg)")
print(f"Resulting speed: {v_rms_zero:.2f} m/s")
print("At absolute zero, there is no thermal motion. The particles are stationary.")
print("Without motion, no particles can cross the door to be trapped on one side.\n")
print("Therefore, a non-zero Temperature is the essential experimental parameter.")
<<<B>>>