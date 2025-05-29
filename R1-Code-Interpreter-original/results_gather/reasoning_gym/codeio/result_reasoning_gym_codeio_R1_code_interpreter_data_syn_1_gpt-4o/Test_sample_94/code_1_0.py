import math

def find_feasible_input():
    # Constants
    g = 9.8  # acceleration due to gravity in m/s^2
    E_steel = 2 * 10**11  # Young's modulus for steel in N/m^2
    E_brass = 0.9 * 10**11  # Young's modulus for brass in N/m^2
    
    # Given elongations
    target_elongation_steel = 4996.0  # in 10^(-6) m
    target_elongation_brass = 8015.2  # in 10^(-6) m
    
    # Iterate over possible values
    for d in range(1, 10):  # diameter in 10^(-4) meters
        for l1 in range(1, 10):  # length of steel wire in meters
            for l2 in range(1, 10):  # length of brass wire in meters
                # Calculate cross-sectional area
                A = math.pi * (d * 10**(-4))**2
                
                for m1 in range(1, 1000):
                    for m2 in range(1, 1000):
                        # Calculate forces
                        F_steel = (m1 + m2) * g
                        F_brass = m2 * g
                        
                        # Calculate elongations
                        e_steel = (F_steel * l1) / (A * E_steel) * 10**6
                        e_brass = (F_brass * l2) / (A * E_brass) * 10**6
                        
                        # Check if they match the target elongations within a small tolerance
                        if abs(e_steel - target_elongation_steel) < 0.1 and abs(e_brass - target_elongation_brass) < 0.1:
                            return {"m1": m1, "m2": m2, "l1": l1, "l2": l2, "d": d}

# Find and print the feasible input
feasible_input = find_feasible_input()
print(feasible_input)