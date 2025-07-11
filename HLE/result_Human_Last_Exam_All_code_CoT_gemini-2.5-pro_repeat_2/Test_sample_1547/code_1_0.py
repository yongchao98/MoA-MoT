import random

def simulate_srw_visits():
    """
    Simulates a simple random walk (SRW) in 3D and counts visits to a specific set A.

    The set A is the plane where the first coordinate is zero.
    For this set A, the condition P_x(tau_A < infinity) = 1 holds for all x.
    The simulation demonstrates that A is not transient by showing that it is
    visited many times during a long walk.
    """
    # Parameters
    dimension = 3
    num_steps = 1_000_000
    
    # Start the walk at the origin (0, 0, 0), which is in the set A.
    position = [0] * dimension
    
    # The set A is defined by the condition: position[0] == 0.
    # We start in A, so the visit count begins at 1.
    visits_to_A = 1
    
    # Available moves: step +1 or -1 along one of the axes
    axes = list(range(dimension))
    directions = [-1, 1]
    
    for step in range(num_steps):
        # Choose a random axis to move along
        axis_to_move = random.choice(axes)
        # Choose a random direction (+1 or -1)
        direction = random.choice(directions)
        
        # Update the position
        position[axis_to_move] += direction
        
        # Check if the new position is in the set A
        if position[0] == 0:
            visits_to_A += 1
            
    print(f"Simulation of SRW in {dimension}D for {num_steps} steps.")
    print(f"The set A is the plane where the first coordinate is 0.")
    print(f"The walk started at the origin, which is in A.")
    print(f"Total number of visits to set A: {visits_to_A}")
    
    # For a transient set, the number of visits should be small and finite.
    # The large number of visits suggests the set is recurrent, not transient.

simulate_srw_visits()
