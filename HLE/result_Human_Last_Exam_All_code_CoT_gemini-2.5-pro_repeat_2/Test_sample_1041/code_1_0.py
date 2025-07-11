import random

def simulate_random_walk():
    """
    Simulates a 2D random walk and calculates the probability of hitting a target set.
    """
    start_pos = (0, 300)
    target_set = {(0, 0), (2, 0)}
    
    # We check against the squared radius to avoid costly square root operations.
    radius_squared = 1000**2
    
    # Number of simulations to run. More simulations yield a more accurate result.
    num_simulations = 100000
    
    success_count = 0
    
    # Possible moves: up, down, left, right
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    
    for i in range(num_simulations):
        # Start each walk from the initial position
        pos = start_pos
        
        # Limit steps to prevent extremely long (though theoretically finite) walks
        for _ in range(1000000): 
            # Check for success condition
            if pos in target_set:
                success_count += 1
                break
            
            # Check for failure condition
            if pos[0]**2 + pos[1]**2 > radius_squared:
                break
                
            # Perform a random move
            move = random.choice(moves)
            pos = (pos[0] + move[0], pos[1] + move[1])

    # Calculate the probability
    probability = success_count / num_simulations
    
    print(f"Starting position: {start_pos}")
    print(f"Target set: {target_set}")
    print(f"Disk radius: 1000")
    print(f"Number of simulations: {num_simulations}")
    print(f"Number of successful walks: {success_count}")
    print(f"Estimated probability: {probability:.3f}")

if __name__ == '__main__':
    simulate_random_walk()
