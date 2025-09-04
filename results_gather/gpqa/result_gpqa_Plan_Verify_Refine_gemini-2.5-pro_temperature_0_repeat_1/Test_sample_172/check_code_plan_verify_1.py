import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    
    Question: If uncertainty in space of electron's location, which is travelling 
    with speed v= 2* 10^8 m/s along x-direction is Δx=0.1 nm . Based on the information 
    estimate the minimum uncertainty in the energy ΔE of electron.
    
    Answer: D) ~10^(-16) J
    """

    # --- Given values and constants ---
    v = 2.0 * 10**8  # Speed of the electron in m/s
    delta_x = 0.1 * 10**-9  # Uncertainty in position in meters (0.1 nm)
    hbar = 1.054571817e-34  # Reduced Planck constant in J*s

    # The selected answer is D, which corresponds to an energy of the order of 10^-16 J.
    expected_order_of_magnitude = -16

    # --- Step 1: Calculate the minimum uncertainty in momentum (Δp_x) ---
    # According to the Heisenberg Uncertainty Principle, Δx * Δp_x >= ħ / 2.
    # The minimum uncertainty is therefore Δp_x = ħ / (2 * Δx).
    try:
        delta_p_x = hbar / (2 * delta_x)
    except ZeroDivisionError:
        return "Constraint failed: The uncertainty in position (Δx) cannot be zero."

    # --- Step 2: Calculate the minimum uncertainty in energy (ΔE) ---
    # The uncertainty in energy can be related to the uncertainty in momentum
    # by ΔE ≈ (dE/dp) * Δp.
    # The derivative dE/dp is the group velocity of the particle, which is its speed v.
    # This holds for both non-relativistic (E = p^2/2m) and relativistic (E^2 = (pc)^2 + (mc^2)^2) cases.
    # Therefore, ΔE ≈ v * Δp_x.
    delta_E = v * delta_p_x

    # --- Step 3: Verify the result ---
    # We check if the calculated ΔE has the same order of magnitude as the answer.
    # The order of magnitude can be found using the base-10 logarithm.
    # For a number like 1.055e-16, log10(1.055e-16) is approx -15.97.
    # The floor of this value gives the exponent for the order of magnitude.
    if delta_E == 0:
        calculated_order_of_magnitude = float('-inf')
    else:
        calculated_order_of_magnitude = math.floor(math.log10(delta_E))

    # Compare the calculated order of magnitude with the one from the answer.
    if calculated_order_of_magnitude == expected_order_of_magnitude:
        return "Correct"
    else:
        return (f"The answer is incorrect. "
                f"The calculated minimum uncertainty in energy is ΔE ≈ {delta_E:.3e} J. "
                f"This corresponds to an order of magnitude of 10^{calculated_order_of_magnitude}, "
                f"but the selected answer D implies an order of magnitude of 10^{expected_order_of_magnitude}.")

# Execute the check
result = check_answer()
print(result)