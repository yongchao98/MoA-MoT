import sys

def calculate_neff_change():
    """
    Illustrates how a new particle decaying into neutrinos affects Neff.

    This is a simplified model to demonstrate the logic, not a full cosmological calculation.
    The energy densities are represented in arbitrary units.
    """

    # In the Standard Model (SM), Neff is ~3.044, corresponding to 3 neutrino species.
    N_eff_SM = 3.044

    # Let's define the energy density of a single standard neutrino species as 1 unit.
    rho_1_nu_species = 1.0
    print(f"Let the energy density of one standard neutrino species be: {rho_1_nu_species:.4f} units.\n")

    # The total energy density in SM neutrinos is proportional to N_eff_SM.
    rho_nu_SM = N_eff_SM * rho_1_nu_species
    print(f"The Standard Model (SM) value is N_eff_SM = {N_eff_SM:.4f}")
    print(f"This corresponds to a total neutrino energy density of rho_nu_SM = {N_eff_SM:.4f} * {rho_1_nu_species:.4f} = {rho_nu_SM:.4f} units.\n")

    # Now, introduce the new particle X. It has "non-negligible abundance",
    # so its energy density, rho_X, is greater than zero.
    # Let's assume its decay injects an extra 0.5 energy units into the neutrinos.
    rho_X_decay_injection = 0.5
    print(f"A new massive particle decays, injecting its energy density into neutrinos.")
    print(f"Energy injected from new particle decay, rho_X = {rho_X_decay_injection:.4f} units.\n")

    # The new total energy density of neutrinos is the sum of the SM part and the injected part.
    rho_nu_new = rho_nu_SM + rho_X_decay_injection
    print("The new total energy density in neutrinos is the sum of the SM energy and the injected energy:")
    print(f"rho_nu_new = rho_nu_SM + rho_X = {rho_nu_SM:.4f} + {rho_X_decay_injection:.4f} = {rho_nu_new:.4f} units.\n")
    
    # The new N_eff is calculated from this new total neutrino energy density.
    N_eff_new = rho_nu_new / rho_1_nu_species
    
    print("The new N_eff is this total energy divided by the energy of a single species:")
    # The final equation requested by the user
    print(f"N_eff_new = rho_nu_new / rho_1_nu_species = {rho_nu_new:.4f} / {rho_1_nu_species:.4f} = {N_eff_new:.4f}\n")
    
    # Compare the new Neff with the SM value.
    if N_eff_new > N_eff_SM:
        print(f"Since {N_eff_new:.4f} > {N_eff_SM:.4f}, the value of N_eff would increase.")
        # We redirect the final answer to stdout for the platform.
        # This part is just for the specific platform this is running on.
        if 'streamlit' in sys.modules:
            # In a real scenario, you'd just return the value or conclusion.
            # Here, we print the magic string for the platform.
            print("<<<increase>>>")
    elif N_eff_new < N_eff_SM:
        print(f"Since {N_eff_new:.4f} < {N_eff_SM:.4f}, the value of N_eff would decrease.")
        if 'streamlit' in sys.modules:
            print("<<<decrease>>>")
    else:
        print("The value of N_eff would not change.")
        if 'streamlit' in sys.modules:
            print("<<<no change>>>")

calculate_neff_change()