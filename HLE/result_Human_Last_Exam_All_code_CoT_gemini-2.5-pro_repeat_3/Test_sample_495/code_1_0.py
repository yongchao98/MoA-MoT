def calculate_energy_loss():
    """
    Calculates the energy loss per centimeter for alpha particles in standard air.
    """
    # --- Given Parameters ---
    initial_energy_E0 = 8.5  # MeV
    total_range_R = 8.3      # cm
    distance_x = 4.0         # cm

    # --- Step 1: Calculate the constant 'k' from the range-energy formula R = k * E^(3/2) ---
    # k = R / E0^(3/2)
    k = total_range_R / (initial_energy_E0**1.5)

    # --- Step 2: Calculate the remaining range and energy at distance x ---
    # Remaining range at distance x
    remaining_range_R_rem = total_range_R - distance_x

    # Energy at distance x, calculated from the remaining range: E(x) = (R_rem / k)^(2/3)
    energy_at_x = (remaining_range_R_rem / k)**(2.0 / 3.0)

    # --- Step 3: Calculate the energy loss per centimeter |dE/dx| ---
    # The formula, derived from the range-energy relation, is |dE/dx| = (2/3) * E(x) / (R - x)
    energy_loss_per_cm = (2.0 / 3.0) * energy_at_x / remaining_range_R_rem

    # --- Print the final equation with all numbers ---
    print(f"The energy loss per centimeter |dE/dx| is calculated as follows:")
    print(f"|dE/dx| = (2/3) * E(at x={distance_x}cm) / (R - x)")
    print(f"|dE/dx| = (2/3) * {energy_at_x:.3f} MeV / ({total_range_R:.1f} cm - {distance_x:.1f} cm)")
    print(f"|dE/dx| = (2/3) * {energy_at_x:.3f} MeV / {remaining_range_R_rem:.1f} cm")
    print(f"Calculated Energy Loss: {energy_loss_per_cm:.3f} MeV/cm")

if __name__ == '__main__':
    calculate_energy_loss()