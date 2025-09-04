import math

def check_astronomy_answer():
    """
    This function checks the correctness of the answer to the given astronomy question.
    It calculates the rest-frame energy of an observed absorption line and identifies
    the corresponding astrophysical source.
    """

    # --- Given Information & Constants ---
    # Distance to the structure in Gigaparsecs (Gpc)
    distance_gpc = 2.1
    # Observed energy of the absorption line in electron volts (eV)
    energy_observed_ev = 3.9e-6
    # The proposed answer to check
    proposed_answer = 'C'  # C) Cold atomic interstellar medium.

    # --- Physical and Cosmological Constants ---
    # Hubble Constant in km/s/Mpc. 70 is a standard approximate value.
    H0_kms_per_mpc = 70.0
    # Speed of light in km/s
    c_kms = 300000.0
    # Planck's constant in eV·s
    h_ev_s = 4.1357e-15
    # Theoretical frequency of the 21-cm hydrogen line in Hz
    freq_21cm_hz = 1420.40575e6

    # --- Step 1: Calculate the cosmological redshift (z) ---
    # Convert distance from Gpc to Mpc for use with H0
    distance_mpc = distance_gpc * 1000
    # Use the Hubble-Lemaître law (v = H0 * D) to find the recessional velocity.
    recessional_velocity_kms = H0_kms_per_mpc * distance_mpc
    # Calculate redshift z = v/c
    redshift = recessional_velocity_kms / c_kms

    # --- Step 2: Calculate the rest-frame energy (E_rest) ---
    # The observed energy is related to the rest-frame energy by E_obs = E_rest / (1 + z)
    # Therefore, E_rest = E_obs * (1 + z)
    energy_rest_calculated_ev = energy_observed_ev * (1 + redshift)

    # --- Step 3: Identify the known transition corresponding to this energy ---
    # The most prominent low-energy transition in the interstellar medium is the 21-cm line
    # of neutral atomic hydrogen (H I). Let's calculate its theoretical energy.
    # Energy = Planck's constant * frequency (E = hf)
    energy_21cm_theoretical_ev = h_ev_s * freq_21cm_hz

    # --- Step 4: Verify the match ---
    # Check if the calculated rest-frame energy matches the 21-cm line energy.
    # A 5% tolerance is used to account for uncertainty in the Hubble Constant.
    if not math.isclose(energy_rest_calculated_ev, energy_21cm_theoretical_ev, rel_tol=0.05):
        return (f"Incorrect. The calculated rest-frame energy ({energy_rest_calculated_ev:.3e} eV) "
                f"does not match the theoretical energy of the 21-cm hydrogen line ({energy_21cm_theoretical_ev:.3e} eV) "
                f"within a 5% tolerance. The redshift was calculated to be z={redshift:.3f} based on H0={H0_kms_per_mpc} km/s/Mpc.")

    # --- Step 5: Match the transition to the correct ISM phase ---
    # The 21-cm line is a hyperfine transition in neutral ATOMIC hydrogen (H I).
    # This rules out options A and B, which are MOLECULAR media.
    # An ABSORPTION line at 21-cm is a classic signature of the COLD Neutral Medium (CNM).
    # The Warm Neutral Medium (WNM) is typically observed via its 21-cm EMISSION.
    # Therefore, "Cold atomic interstellar medium" is the most precise and correct description.

    if proposed_answer == 'C':
        return "Correct"
    else:
        return (f"Incorrect. The calculated rest-frame energy corresponds to the 21-cm line of neutral atomic hydrogen. "
                f"An absorption feature at this wavelength is characteristic of the 'Cold atomic interstellar medium' (Option C), "
                f"not the proposed answer '{proposed_answer}'.")

# Execute the check
result = check_astronomy_answer()
print(result)