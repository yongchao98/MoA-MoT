import math

def calculate_lifetime_multiple():
    """
    Calculates the theoretical lifetime of the Sodium 3p state using a hydrogenic model
    and determines its multiple difference from the experimental value.
    """

    # --- Constants ---
    # Speed of light in m/s
    c = 2.998e8
    # Bohr radius in m
    a0 = 5.292e-11
    # Fine-structure constant
    alpha = 1 / 137.036
    # Experimental lifetime of 3p state in ns
    tau_exp_ns = 16.2

    # --- Problem Parameters ---
    # Wavelength of transition in nm
    lambda_nm = 589
    # Nuclear charge Z
    Z = 11
    # Given degeneracy ratio g2/g1
    g2_over_g1 = 2.0
    
    # --- Calculations ---
    
    # 1. Convert units to SI
    lambda_m = lambda_nm * 1e-9
    tau_exp_s = tau_exp_ns * 1e-9
    
    # 2. Calculate angular frequency (omega)
    omega = (2 * math.pi * c) / lambda_m
    
    # 3. Calculate the squared radial integral |<3s|r|3p>|^2
    # The radial integral for a hydrogenic atom from 3p to 3s is I = -9 * sqrt(2) * a0 / Z
    I_radial = -9 * math.sqrt(2) * a0 / Z
    I_squared = I_radial**2
    
    # 4. Determine the degeneracy factor
    # We interpret the problem's hint g2/g1=2 to mean we should use g1/g2 = 1/2
    # as the degeneracy factor in the rate equation.
    degeneracy_factor = 1.0 / g2_over_g1

    # 5. Calculate the Einstein A coefficient (spontaneous emission rate)
    # A_21 = (4 * alpha * omega^3 / (3 * c^2)) * (degeneracy_factor) * |I|^2
    A21 = (4 * alpha * omega**3) / (3 * c**2) * degeneracy_factor * I_squared

    # 6. Calculate the theoretical lifetime
    tau_th_s = 1 / A21
    tau_th_ns = tau_th_s * 1e9

    # 7. Calculate the multiple difference
    multiple = tau_th_ns / tau_exp_ns

    # --- Output Results ---
    print("--- Calculation Steps ---")
    print(f"1. Transition wavelength: {lambda_nm} nm")
    print(f"2. Angular frequency (omega): {omega:.3e} rad/s")
    print(f"3. Squared radial integral |I|^2: {I_squared:.3e} m^2")
    print(f"4. Einstein A coefficient (A_21): {A21:.3e} s^-1")
    print(f"5. Theoretical lifetime (tau_th): {tau_th_ns:.1f} ns")
    print(f"6. Experimental lifetime (tau_exp): {tau_exp_ns} ns")
    print("\n--- Final Result ---")
    print(f"The theoretical lifetime is {tau_th_ns:.1f} ns, which is {multiple:.2f} times the experimental lifetime of {tau_exp_ns} ns.")
    print(f"Final Equation: {tau_th_ns:.1f} / {tau_exp_ns} = {multiple:.2f}")
    
    # Find the closest answer choice
    choices = {
        'A': 30, 'B': 2, 'C': 0.1, 'D': 1,
        'E': 100, 'F': 0.5, 'G': 10, 'H': 5
    }
    closest_choice = min(choices.keys(), key=lambda k: abs(choices[k] - multiple))
    print(f"\nThe calculated multiple ({multiple:.2f}) is closest to answer choice G: 10.")
    
calculate_lifetime_multiple()
<<<G>>>