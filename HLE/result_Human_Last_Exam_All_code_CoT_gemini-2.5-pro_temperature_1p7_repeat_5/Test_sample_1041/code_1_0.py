import math

def calculate_hitting_probability():
    """
    Calculates the probability that a 2D random walk visits a target set before exiting a disk.
    """
    # Problem parameters
    s0 = (0, 300)
    w1 = (0, 0)
    w2 = (2, 0)
    R = 1000
    
    # Effective radius of a lattice point
    rho = 1

    print("Step 1: Calculate individual hitting probabilities P1 and P2.")

    # Distance from s0 to w1 and w2
    dist_s0_w1 = math.sqrt((s0[0] - w1[0])**2 + (s0[1] - w1[1])**2)
    dist_s0_w2 = math.sqrt((s0[0] - w2[0])**2 + (s0[1] - w2[1])**2)

    # Calculate P1 and P2 using the formula P = ln(R/dist) / ln(R/rho)
    p1 = math.log(R / dist_s0_w1) / math.log(R / rho)
    p2 = math.log(R / dist_s0_w2) / math.log(R / rho)
    
    print(f"Start S0 = {s0}, Target w1 = {w1}, Target w2 = {w2}, Radius R = {R}")
    print(f"Distance from S0 to w1: {dist_s0_w1:.4f}")
    print(f"Distance from S0 to w2: {dist_s0_w2:.4f}")
    print(f"Probability P1 (hitting w1 first): ln({R}/{dist_s0_w1:.4f}) / ln({R}) = {p1:.5f}")
    print(f"Probability P2 (hitting w2 first): ln({R}/{dist_s0_w2:.4f}) / ln({R}) = {p2:.5f}")
    
    print("\nStep 2: Calculate the transfer probability q.")
    
    # Neighbors of w1
    neighbors_w1 = [(w1[0] + dx, w1[1] + dy) for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]]
    
    q = 0
    print(f"Calculating q = P(w1 -> w2 before exit) by averaging from neighbors of w1 = {w1}:")
    for v in neighbors_w1:
        dist_v_w2 = math.sqrt((v[0] - w2[0])**2 + (v[1] - w2[1])**2)
        # Hitting probability from a neighbor v to w2
        p_v_w2 = math.log(R / dist_v_w2) / math.log(R / rho)
        q += (1/4) * p_v_w2
        print(f"  Neighbor v = {v}, dist(v, w2) = {dist_v_w2:.4f}, P(v -> w2) = {p_v_w2:.5f}")
        
    print(f"The total transfer probability q = {q:.5f}")
    
    print("\nStep 3: Calculate the final probability of hitting the set A = {w1, w2}.")

    # Final probability P(A) = (P1 + P2) / (1 + q)
    final_prob = (p1 + p2) / (1 + q)

    print(f"Final Probability P(A) = (P1 + P2) / (1 + q)")
    print(f"P(A) = ({p1:.5f} + {p2:.5f}) / (1 + {q:.5f})")
    print(f"P(A) = {final_prob:.5f}")
    
    print(f"\nThe probability is approximately {final_prob:.3f} when rounded to three significant digits.")
    
    
calculate_hitting_probability()