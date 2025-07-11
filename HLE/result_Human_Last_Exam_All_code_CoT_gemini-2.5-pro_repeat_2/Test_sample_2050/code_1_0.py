import math

def calculate_ground_state_energy():
    """
    Calculates the ground state energy for the given system.

    The ground state energy is E_ground = c_omega * omega + c_J * J.
    This function calculates the coefficients c_omega and c_J.
    """
    N = 7  # Number of cavities
    n_photons = 4  # Number of photons

    # The coefficient of omega is just the number of photons.
    c_omega = n_photons

    # The U -> infinity limit means we have hard-core bosons, equivalent to
    # non-interacting fermions. We fill the lowest 4 single-particle energy levels.
    # The single-particle energies are E_m = -2J*cos(2*pi*m/N).
    # We need to find the 4 quantum numbers 'm' that give the lowest energies.
    # The energies sorted from lowest to highest correspond to m = 0, (+1, -1), (+2, -2), ...
    # So we occupy the states with m = 0, 1, -1, 2.
    m_values = [0, 1, -1, 2]

    # The energy from the J term is the sum of the single-particle energies
    # E_J = sum(-2J * cos(2*pi*m/N)) for m in m_values
    # c_J = -2 * sum(cos(2*pi*m/N))
    
    cos_terms = [math.cos(2 * math.pi * m / N) for m in m_values]
    sum_of_cos = sum(cos_terms)
    c_J = -2 * sum_of_cos

    # --- Output the results step-by-step ---
    print("The ground state energy E_ground is of the form: c_omega * omega + c_J * J")
    print("-" * 50)
    
    # Output omega coefficient
    print(f"The number of photons is {n_photons}.")
    print(f"The coefficient of omega is c_omega = {c_omega}")
    print("-" * 50)

    # Output J coefficient details
    print("To find the coefficient of J, we sum the energies of the 4 lowest states.")
    print(f"These states correspond to momentum quantum numbers m = {m_values}.")
    
    # Explicitly show the numbers in the calculation
    cos_0 = math.cos(2 * math.pi * 0 / N)
    cos_2pi_7 = math.cos(2 * math.pi * 1 / N)
    cos_4pi_7 = math.cos(2 * math.pi * 2 / N)
    
    print("The coefficient c_J is given by the expression:")
    print(f"c_J = -2 * (cos(2*pi*0/7) + cos(2*pi*1/7) + cos(2*pi*(-1)/7) + cos(2*pi*2/7))")
    print(f"c_J = -2 * (1 + 2 * cos(2*pi/7) + cos(4*pi/7))")
    print("\nNumerical values for the cosine terms:")
    print(f"cos(0) = {cos_0:.8f}")
    print(f"cos(2*pi/7) = {cos_2pi_7:.8f}")
    print(f"cos(4*pi/7) = {cos_4pi_7:.8f}")
    
    print(f"\nSum of cosine terms = {cos_0:.4f} + 2 * {cos_2pi_7:.4f} + {cos_4pi_7:.4f} = {sum_of_cos:.8f}")
    print(f"The coefficient of J is c_J = -2 * {sum_of_cos:.8f} = {c_J:.8f}")
    print("-" * 50)

    # Final equation
    print("The final equation for the ground state energy is:")
    print(f"E_ground = {c_omega} * omega + ({c_J:.8f}) * J")

calculate_ground_state_energy()