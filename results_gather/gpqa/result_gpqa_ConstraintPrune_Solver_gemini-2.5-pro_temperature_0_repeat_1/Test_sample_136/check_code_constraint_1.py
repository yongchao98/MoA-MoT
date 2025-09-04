import math

def check_synchrocyclotron_revolutions():
    """
    This function checks the correctness of the given answer for the number of revolutions
    a proton makes in a synchrocyclotron.
    """
    # --- Given Data ---
    # Final kinetic energy in Mega-electron-volts (MeV)
    T_final_MeV = 950.0
    # Peak accelerating voltage in kilo-volts (kV)
    U0_kV = 190.0
    # Synchronous phase in radians
    phi0_rad = math.pi / 4
    # The answer to be checked, from option A
    given_answer = 3536

    # --- Physics Calculation ---

    # 1. Convert units to a consistent system (eV and V)
    T_final_eV = T_final_MeV * 1e6  # Final energy in eV
    U0_V = U0_kV * 1e3            # Peak voltage in V

    # 2. Calculate the energy gained in a single acceleration event.
    # The energy gain in eV is numerically equal to the potential difference in Volts
    # at the moment of acceleration.
    # Voltage at acceleration = U0 * cos(phi0)
    energy_gain_per_acceleration_eV = U0_V * math.cos(phi0_rad)

    # 3. Calculate the energy gained per revolution.
    # The problem states there are two accelerations per revolution.
    energy_gain_per_revolution_eV = 2 * energy_gain_per_acceleration_eV

    # 4. Calculate the total number of revolutions.
    # This is the total energy gained divided by the energy gain per revolution.
    # We assume the initial kinetic energy is negligible.
    if energy_gain_per_revolution_eV <= 0:
        return "Error: Energy gain per revolution is non-positive. Cannot calculate revolutions."
    
    calculated_revolutions = T_final_eV / energy_gain_per_revolution_eV

    # --- Verification ---

    # The number of revolutions should be an integer. We check if the given answer
    # is the rounded value of our calculation.
    if round(calculated_revolutions) == given_answer:
        return "Correct"
    else:
        reason = (
            f"The calculated number of revolutions is incorrect.\n"
            f"Calculation Steps:\n"
            f"1. Energy gain per acceleration = U0 * cos(Φ0) = {U0_V:.0f} V * cos(π/4) ≈ {energy_gain_per_acceleration_eV:.2f} eV\n"
            f"2. Energy gain per revolution = 2 * {energy_gain_per_acceleration_eV:.2f} eV ≈ {energy_gain_per_revolution_eV:.2f} eV\n"
            f"3. Total revolutions = Total Energy / Energy per Revolution = {T_final_eV:.0f} eV / {energy_gain_per_revolution_eV:.2f} eV ≈ {calculated_revolutions:.2f}\n"
            f"The calculated value rounds to {round(calculated_revolutions)}, but the given answer is {given_answer}."
        )
        return reason

# The final output will be the result of this function call.
print(check_synchrocyclotron_revolutions())