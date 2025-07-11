def calculate_new_neff(n_eff_sm, rho_x_over_rho_gamma):
    """
    Calculates the change in the effective number of neutrino species (N_eff)
    due to an energy injection from a new particle decaying into neutrinos.

    Args:
        n_eff_sm (float): The Standard Model value of N_eff.
        rho_x_over_rho_gamma (float): The energy density of the new particle (X)
                                     as a fraction of the photon energy density
                                     just before it decays.
    """
    # N_eff parameterizes the radiation energy density:
    # rho_rad = rho_gamma * (1 + N_eff * (7/8) * (T_nu/T_gamma)^4)
    # The term (T_nu/T_gamma)^4 is (4/11)^(4/3) after e+e- annihilation.
    # The change in N_eff is related to the injected neutrino energy density (delta_rho_nu)
    # delta_rho_nu = delta_N_eff * rho_gamma * (7/8) * (4/11)^(4/3)
    # The injected energy comes from the new particle, so delta_rho_nu = rho_x.
    # Rearranging for delta_N_eff:
    # delta_N_eff = (rho_x / rho_gamma) / ((7/8) * (4/11)^(4/3))

    print(f"Starting with the Standard Model value N_eff_SM = {n_eff_sm:.4f}\n")

    # Define the ratio constant from the energy density formula
    # k = (7/8) * (T_nu / T_gamma)^4
    k = (7.0 / 8.0) * (4.0 / 11.0)**(4.0 / 3.0)

    # Let's assume the energy density of the new particle X is 1% of the photon energy density
    print(f"Assuming the injected energy density from the new particle is {rho_x_over_rho_gamma*100:.1f}% of the photon energy density.")
    print(f"This gives a ratio (rho_X / rho_gamma) = {rho_x_over_rho_gamma:.4f}\n")
    
    # Calculate the change in N_eff
    delta_n_eff = rho_x_over_rho_gamma / k
    
    # Calculate the new N_eff
    n_eff_new = n_eff_sm + delta_n_eff

    print("The final N_eff is calculated as: N_eff_new = N_eff_SM + delta_N_eff")
    print("The increase, delta_N_eff, is calculated from the injected energy.")
    print(f"delta_N_eff = ({rho_x_over_rho_gamma:.4f}) / ((7/8) * (4/11)^(4/3))")
    print(f"delta_N_eff = ({rho_x_over_rho_gamma:.4f}) / ({k:.4f}) = {delta_n_eff:.4f}\n")
    
    print("The final result is:")
    print(f"{n_eff_new:.4f} = {n_eff_sm:.4f} + {delta_n_eff:.4f}")

    if n_eff_new > n_eff_sm:
        print("\nAs shown by the calculation, N_eff has increased.")
    else:
        # This case should not be reached with a positive energy injection.
        print("\nN_eff has not increased.")

if __name__ == '__main__':
    # The standard model value for N_eff, including small corrections.
    N_eff_Standard_Model = 3.046
    
    # Let's assume the decaying particle's energy density is a small fraction
    # (e.g., 1%) of the photon energy density at the time of decay.
    rho_X_relative_to_rho_gamma = 0.01

    calculate_new_neff(N_eff_Standard_Model, rho_X_relative_to_rho_gamma)