import sys

def calculate_efficiency():
    """
    This function models and compares the quantum efficiency of a non-radical
    and a radical-based OLED emitter, highlighting the effect of self-quenching.
    """
    # --- Shared Parameters ---
    # Radiative decay rate (how fast an excited state emits a photon).
    k_r = 1.0e7  # in s^-1
    # Intrinsic non-radiative decay rate (e.g., from vibrations).
    k_nr = 1.0e6  # in s^-1
    # Assumed concentration of the emitter in the emissive layer.
    concentration = 0.1  # in M (moles/liter)
    # Outcoupling efficiency (fraction of photons that escape the device).
    eta_outcoupling = 0.25

    # --- Specific Parameters ---
    # For a typical non-radical material, bimolecular self-quenching is negligible.
    k_q_non_radical = 0.0
    # For a radical, the quenching rate constant (k_q) is significant because
    # an excited radical can be quenched by a ground-state radical.
    k_q_radical = 5.0e8  # in M^-1 s^-1

    print("--- OLED Efficiency Model ---")
    print("This model demonstrates the impact of self-quenching on radical emitters.")
    print("The key equation is: EQE = eta_outcoupling * k_r / (k_r + k_nr + k_q * [Q])\n")

    # --- Calculation for Non-Radical Emitter ---
    denominator_non_radical = k_r + k_nr + k_q_non_radical * concentration
    iqe_non_radical = k_r / denominator_non_radical
    eqe_non_radical = iqe_non_radical * eta_outcoupling

    print("--- Case 1: Non-Radical Emitter ---")
    print("The final equation for the non-radical emitter's EQE is:")
    # Using 'g' format specifier for clean scientific notation
    print(f"EQE = {eta_outcoupling} * {k_r:.1g} / ({k_r:.1g} + {k_nr:.1g} + {k_q_non_radical:.1g} * {concentration})")
    print(f"Resulting EQE: {eqe_non_radical:.2%}\n")


    # --- Calculation for Radical Emitter ---
    denominator_radical = k_r + k_nr + k_q_radical * concentration
    iqe_radical = k_r / denominator_radical
    eqe_radical = iqe_radical * eta_outcoupling

    print("--- Case 2: Radical Emitter (with self-quenching) ---")
    print("The final equation for the radical emitter's EQE is:")
    print(f"EQE = {eta_outcoupling} * {k_r:.1g} / ({k_r:.1g} + {k_nr:.1g} + {k_q_radical:.1g} * {concentration})")
    print(f"Resulting EQE: {eqe_radical:.2%}\n")

    print("Conclusion: The radical-radical quenching pathway (k_q) significantly reduces the number of emitted photons, leading to a much lower EQE.")

if __name__ == '__main__':
    calculate_efficiency()