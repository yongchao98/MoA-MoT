import math

def calculate_u_r_sequence():
    """
    Calculates the sequence u_r(n) for n from 3 to 12.
    
    The term u_r(n) denotes the minimal order of the Picard-Fuchs differential
    equation for the real-analytic periods on alpha in (0, infinity)
    for the Hamiltonian V(q) = 1/2 * (q^2 - q^n).
    
    The formula used is u_r(n) = floor(n/2) - 1, based on research by
    Konovalov, Ternovsky, and Vashchenko.
    """
    results = []
    
    print("This program calculates the sequence {u_r(3), u_r(4), ..., u_r(12)}.")
    print("The formula is u_r(n) = floor(n/2) - 1.")
    print("-" * 50)

    for n in range(3, 13):
        # Calculate the value using the formula
        n_over_2 = n / 2
        floor_val = math.floor(n_over_2)
        u_val = floor_val - 1
        results.append(u_val)
        
        # Print the detailed calculation for each n
        print(f"For n = {n}:")
        print(f"u_r({n}) = floor({n}/2) - 1 = floor({n_over_2}) - 1 = {floor_val} - 1 = {u_val}")
        print("-" * 50)
        
    print("The complete sequence {u_r(3), u_r(4), ..., u_r(12)} is:")
    print(results)

# Execute the function
calculate_u_r_sequence()