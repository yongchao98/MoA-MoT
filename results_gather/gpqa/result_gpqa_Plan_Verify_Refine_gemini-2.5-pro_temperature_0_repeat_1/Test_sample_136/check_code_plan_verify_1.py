import math

def check_answer():
    """
    This function verifies the answer to the synchrocyclotron problem.
    It calculates the number of revolutions a proton makes to reach a
    final kinetic energy based on the given parameters.
    """
    # --- Define problem constants and the given answer ---
    # Final kinetic energy in MeV
    T_final_MeV = 950.0
    # Peak accelerating voltage in kV
    U0_kV = 190.0
    # Synchronous phase in radians
    phi0 = math.pi / 4
    # The answer from the LLM to be checked (Option C)
    llm_answer_value = 3536

    # --- Convert to consistent base units (eV and V) ---
    # Final kinetic energy in electron-volts (eV)
    T_final_eV = T_final_MeV * 1e6
    # Peak voltage in Volts (V)
    U0_V = U0_kV * 1e3

    # --- Perform the physics calculation ---
    # The energy gained by a proton (charge q=1e) when crossing a potential
    # difference U (in Volts) is ΔT = U (in eV).
    # In a synchrocyclotron, the particle is accelerated twice per revolution.
    # The voltage at the moment of crossing is assumed to be at the synchronous
    # phase: U_sync = U0 * cos(phi0).
    # Therefore, the energy gained per revolution is ΔT_rev = 2 * U_sync.
    
    try:
        # Calculate the energy gained per revolution in eV
        energy_gain_per_revolution_eV = 2 * U0_V * math.cos(phi0)

        if energy_gain_per_revolution_eV <= 0:
            return "Calculation error: Energy gained per revolution must be positive."

        # The total number of revolutions (N) is the total energy gained divided
        # by the energy gained per revolution. Initial kinetic energy is negligible.
        # N = T_final / ΔT_rev
        calculated_revolutions = T_final_eV / energy_gain_per_revolution_eV

        # The number of revolutions should be an integer. We round the result to the nearest whole number.
        rounded_calculated_revolutions = round(calculated_revolutions)

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Compare the calculated result with the LLM's answer ---
    if rounded_calculated_revolutions == llm_answer_value:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        reason = (
            f"The provided answer {llm_answer_value} is incorrect.\n"
            f"The calculation shows a different result.\n"
            f"1. The total kinetic energy to be gained is T_final = {T_final_MeV} MeV = {T_final_eV:.2e} eV.\n"
            f"2. The energy gained per revolution is ΔT_rev = 2 * U₀ * cos(Φ₀).\n"
            f"   With U₀ = {U0_kV} kV and Φ₀ = π/4, ΔT_rev = 2 * {U0_V:.0f} V * cos(π/4) ≈ {energy_gain_per_revolution_eV:.2f} eV.\n"
            f"3. The total number of revolutions is N = T_final / ΔT_rev.\n"
            f"   N = {T_final_eV:.2e} eV / {energy_gain_per_revolution_eV:.2f} eV ≈ {calculated_revolutions:.4f}.\n"
            f"4. Rounding this to the nearest integer gives {rounded_calculated_revolutions}.\n"
            f"5. This calculated value ({rounded_calculated_revolutions}) does not match the provided answer ({llm_answer_value})."
        )
        return reason

# To run the check, you would call the function and print its output:
# result = check_answer()
# print(result)