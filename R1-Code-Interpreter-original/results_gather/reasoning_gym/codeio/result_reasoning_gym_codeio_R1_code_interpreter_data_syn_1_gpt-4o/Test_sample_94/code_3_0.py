import math

def calculate_masses():
    # Constants
    g = 9.8  # acceleration due to gravity in m/s^2
    E_steel = 2 * 10**11  # Young's modulus for steel in N/m^2
    E_brass = 0.9 * 10**11  # Young's modulus for brass in N/m^2
    
    # Given elongations
    target_elongation_steel = 4996.0  # in 10^(-6) m
    target_elongation_brass = 8015.2  # in 10^(-6) m
    
    # Assume reasonable values for l1, l2, and d
    d = 5  # diameter in 10^(-4) meters
    l1 = 2  # length of steel wire in meters
    l2 = 3  # length of brass wire in meters
    
    # Calculate cross-sectional area
    A = math.pi * (d * 10**(-4))**2
    
    # Calculate the required forces based on the target elongations
    F_steel = (target_elongation_steel / 10**6) * (A * E_steel) / l1
    F_brass = (target_elongation_brass / 10**6) * (A * E_brass) / l2
    
    # Calculate the masses
    m1_plus_m2 = F_steel / g
    m2 = F_brass / g
    m1 = m1_plus_m2 - m2
    
    return {"m1": round(m1), "m2": round(m2), "l1": l1, "l2": l2, "d": d}

# Calculate and print the masses
feasible_input = calculate_masses()
print(feasible_input)