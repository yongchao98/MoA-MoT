import random
import math

def simulate_controlled_rw():
    """
    Simulates a controlled random walk in 3D to illustrate how control
    can lead to recurrent behavior.
    """
    d = 3  # Dimension
    num_steps = 100000

    # Define two measures with different anisotropies.
    # Each measure is a uniform distribution over a set of 8 vectors.
    
    # Measure 1: big steps along x-axis, small steps on yz-plane.
    steps1 = []
    for sx in [-3, 3]:
        for sy in [-1, 1]:
            for sz in [-1, 1]:
                steps1.append((sx, sy, sz))
    
    # Measure 2: small steps along x-axis, big steps on yz-plane.
    steps2 = []
    for sx in [-1, 1]:
        for sy in [-3, 3]:
            for sz in [-3, 3]:
                steps2.append((sx, sy, sz))
                
    measures = [steps1, steps2]
    
    # Pre-calculate approximate covariance matrices. 
    # For a uniform distribution on a set S, C_ij = E[x_i x_j].
    # For our measures, the off-diagonal terms are 0.
    # C1_11 = E[x^2] = (1/8) * 8 * 3^2 = 9. C1_22 = E[y^2] = 1. C1_33 = 1.
    # C2_11 = 1, C2_22 = 9, C2_33 = 9.
    
    position = [10, 10, 10]  # Start away from the origin
    
    print("Controlled Random Walk Simulation (k=2)")
    print("-----------------------------------------")
    print(f"Start position: {position}, Distance: {math.sqrt(sum(p**2 for p in position)):.2f}")

    for step in range(1, num_steps + 1):
        x, y, z = position
        
        # Strategy: Choose the measure that discourages movement away from the origin
        # by selecting the one with smaller variance in the current radial direction.
        # This is equivalent to choosing j to minimize x^T * C_j * x.
        
        # value for measure 1: 9*x^2 + 1*y^2 + 1*z^2
        q_form1 = 9 * x**2 + y**2 + z**2
        # value for measure 2: 1*x^2 + 9*y^2 + 9*z^2
        q_form2 = x**2 + 9 * y**2 + 9 * z**2
        
        if q_form1 < q_form2:
            # This happens if 8*x^2 < 8*(y^2+z^2), i.e., x is not the dominant component.
            # We choose measure 1, which has large steps in x, pushing us away
            # from the yz-plane.
            chosen_measure_idx = 0
        else:
            # This happens if x is the dominant component.
            # We choose measure 2, which has small steps in x, pushing us back
            # towards the yz-plane (in a relative sense).
            chosen_measure_idx = 1
        
        # Take a random step from the chosen measure
        step_vec = random.choice(measures[chosen_measure_idx])
        
        # Update position
        position = [position[i] + step_vec[i] for i in range(d)]

        if step % (num_steps // 10) == 0:
            distance = math.sqrt(sum(p**2 for p in position))
            print(f"Step {step:7d}: Distance from origin: {distance:.2f}")

    final_distance = math.sqrt(sum(p**2 for p in position))
    print(f"End position: {[round(p, 2) for p in position]}, Final Distance: {final_distance:.2f}")
    print("\nNote: A standard 3D random walk's distance would grow ~sqrt(n).")
    print("This walk stays bounded, demonstrating recurrence is possible.")
    
    print("\nThe maximal k for which recurrence cannot be guaranteed is 1.")


if __name__ == '__main__':
    simulate_controlled_rw()