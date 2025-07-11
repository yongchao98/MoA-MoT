import math

def solve_random_walk_probability():
    """
    Calculates the probability for a 2D random walk starting at (0, 300) to visit
    the set {(0,0), (2,0)} before leaving the disk of radius 1000.
    """
    
    # Parameters of the problem
    R = 1000.0  # Radius of the disk
    z0 = (0, 300)  # Starting point
    w1 = (0, 0)    # First target point
    w2 = (2, 0)    # Second target point
    
    print("Step 1: Define the parameters.")
    print(f"Radius of the disk, R = {R}")
    print(f"Starting point, z0 = {z0}")
    print(f"Target set, A = {{ {w1}, {w2} }}")
    print("-" * 30)

    # Use the approximate Green's function for a large disk:
    # G_R(z, w) = log(R) - log|z-w| for z != w
    # G_R(w, w) = log(R)
    
    # We set up a system of linear equations to find the "charges" e1 and e2
    # on w1 and w2 such that the potential on the target set is 1.
    # [ G(w1,w1) G(w1,w2) ] [e1] = [1]
    # [ G(w2,w1) G(w2,w2) ] [e2]   [1]

    log_R = math.log(R)
    
    # Off-diagonal term G(w1, w2) = G(w2, w1)
    dist_w1_w2 = math.sqrt((w1[0] - w2[0])**2 + (w1[1] - w2[1])**2)
    G_w1_w2 = log_R - math.log(dist_w1_w2)
    
    # Diagonal terms are approximately log(R)
    G_w1_w1 = log_R
    G_w2_w2 = log_R
    
    print("Step 2: Set up the linear system for the equilibrium charges e1, e2.")
    print(f"G(w1,w1) = log(R) = {G_w1_w1:.4f}")
    print(f"G(w2,w2) = log(R) = {G_w2_w2:.4f}")
    print(f"G(w1,w2) = log(R) - log(|w1-w2|) = log(1000) - log({dist_w1_w2}) = {G_w1_w2:.4f}")
    
    # The system is:
    # e1 * log_R + e2 * (log_R - log(2)) = 1
    # e1 * (log_R - log(2)) + e2 * log_R = 1
    # By symmetry of the equations, e1 = e2. Let's call it 'e'.
    # e * (log_R + log_R - log(2)) = 1
    
    e_denominator = 2 * log_R - math.log(dist_w1_w2)
    e = 1.0 / e_denominator
    e1, e2 = e, e
    
    print(f"\nSolving the system gives e1 = e2 = e.")
    print(f"e = 1 / (2*log(R) - log(2)) = 1 / {e_denominator:.4f} = {e:.4f}")
    print("-" * 30)
    
    # The probability is the potential at z0:
    # h(z0) = e1 * G(z0, w1) + e2 * G(z0, w2)
    
    dist_z0_w1 = math.sqrt((z0[0] - w1[0])**2 + (z0[1] - w1[1])**2)
    dist_z0_w2 = math.sqrt((z0[0] - w2[0])**2 + (z0[1] - w2[1])**2)
    
    G_z0_w1 = log_R - math.log(dist_z0_w1)
    G_z0_w2 = log_R - math.log(dist_z0_w2)

    print("Step 3: Calculate the probability (potential) at the starting point z0.")
    print(f"Distance from z0 to w1: |z0-w1| = {dist_z0_w1:.4f}")
    print(f"Distance from z0 to w2: |z0-w2| = {dist_z0_w2:.4f}")

    # Final probability calculation
    numerator = (G_z0_w1 + G_z0_w2)
    probability = e * numerator
    
    print(f"\nThe probability is given by the formula:")
    print(f"P = e * (G(z0,w1) + G(z0,w2))")
    print(f"P = e * ( (log(R) - log(|z0-w1|)) + (log(R) - log(|z0-w2|)) )")
    print(f"P = e * ( 2*log(R) - log(|z0-w1|) - log(|z0-w2|) )")
    
    final_numerator = 2*log_R - math.log(dist_z0_w1) - math.log(dist_z0_w2)
    final_denominator = e_denominator
    
    print(f"\nFinal calculation:")
    print(f"Numerator = 2*log({R}) - log({dist_z0_w1:.4f}) - log({dist_z0_w2:.4f}) = {final_numerator:.4f}")
    print(f"Denominator = 2*log({R}) - log({dist_w1_w2}) = {final_denominator:.4f}")
    print(f"Probability = Numerator / Denominator = {final_numerator:.4f} / {final_denominator:.4f} = {probability:.4f}")
    
    print("\n" + "="*50)
    print(f"The probability that the random walk visits the set before leaving the disk is approximately: {probability:.3f}")
    print("="*50)

solve_random_walk_probability()