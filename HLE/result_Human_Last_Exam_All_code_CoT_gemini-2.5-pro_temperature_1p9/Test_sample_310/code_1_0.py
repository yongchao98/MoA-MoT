import sys

def calculate_proton_beam_properties(total_available_energy_MeV, oam_mode):
    """
    Models the effect of Orbital Angular Momentum (OAM) on a proton beam.

    Args:
        total_available_energy_MeV (float): The maximum possible proton energy without OAM (l=0).
        oam_mode (int): The OAM mode number (l). l=0 for a standard beam.
    """

    # --- Model Parameters ---
    # Assume 10% of energy is partitioned into rotation for each OAM mode 'l'.
    energy_partition_factor = 0.10
    # Assume a base beam dispersion of 5 mrad, which increases by 3 mrad per OAM mode.
    base_dispersion_mrad = 5.0
    dispersion_factor_mrad = 3.0

    print("--- Proton Beam Property Calculation ---")
    print(f"Scenario: Laser with OAM mode l = {oam_mode}")
    print(f"Total energy available from laser pulse for proton acceleration: {total_available_energy_MeV} MeV\n")

    # 1. Calculate Energy Partitioning
    print("Step 1: Calculating Proton Energy")
    print("Model Equation: Final_Proton_Energy = Total_Energy - Rotational_Energy")

    rotational_energy_MeV = total_available_energy_MeV * energy_partition_factor * oam_mode
    final_proton_energy_MeV = total_available_energy_MeV - rotational_energy_MeV

    # Outputting the numbers in the final equation
    print("Result:")
    print(f"{final_proton_energy_MeV:.2f} MeV = {total_available_energy_MeV:.2f} MeV - {rotational_energy_MeV:.2f} MeV")
    if oam_mode > 0:
        print("Conclusion: Final proton energy has DECREASED.\n")
    else:
        print("Conclusion: No energy partitioned to rotation.\n")


    # 2. Calculate Beam Dispersion
    print("Step 2: Calculating Beam Dispersion")
    print("Model Equation: Final_Dispersion = Base_Dispersion + OAM_Induced_Dispersion")

    oam_induced_dispersion = dispersion_factor_mrad * oam_mode
    final_dispersion_mrad = base_dispersion_mrad + oam_induced_dispersion

    # Outputting the numbers in the final equation
    print("Result:")
    print(f"{final_dispersion_mrad:.2f} mrad = {base_dispersion_mrad:.2f} mrad + {oam_induced_dispersion:.2f} mrad")

    if oam_mode > 0:
        print("Conclusion: Beam DISPERSION has INCREASED.\n")
    else:
        print("Conclusion: Beam has its natural base dispersion.\n")

# --- Main execution ---
if __name__ == "__main__":
    # Case 1: Standard laser beam without OAM (l=0)
    # calculate_proton_beam_properties(total_available_energy_MeV=50.0, oam_mode=0)

    # Case 2: Laser beam with OAM (l=2)
    # This is the case described in the question.
    calculate_proton_beam_properties(total_available_energy_MeV=50.0, oam_mode=2)

    # Final summary based on the model
    print("Summary of effects when imbuing photons with OAM:")
    print("- Energy is partitioned into rotation, reducing the forward kinetic energy of protons.")
    print("- The imparted angular momentum causes the proton beam to spread out, or disperse.")
    sys.stdout.flush()
