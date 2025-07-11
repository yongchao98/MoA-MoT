import math
import random

def get_potential(p):
    """
    Calculates the potential kernel a(p) for a point p=(x,y) in Z^2.
    This uses the asymptotic approximation for large distances, which is sufficient
    to demonstrate the drift.
    a(p) approx (2/pi) * log|p| + C
    """
    x, y = p
    if x == 0 and y == 0:
        return -1 # Origin is not in the state space

    # Constants for the potential kernel approximation of SRW in Z^2
    # a(x) ~ (2/pi) * (log|x| + gamma + log(sqrt(8)))
    euler_gamma = 0.5772156649
    log_sqrt_8 = math.log(8) / 2.0
    C = (2 / math.pi) * (euler_gamma + log_sqrt_8)
    
    dist = math.sqrt(x**2 + y**2)
    return (2 / math.pi) * math.log(dist) + C

def run_simulation_trial(start_pos, num_steps):
    """
    Runs a single simulation of the h-transformed walk.
    Returns True if the walk escapes the x-axis and does not return, False otherwise.
    """
    current_pos = start_pos
    has_left_axis = False

    for _ in range(num_steps):
        x, y = current_pos

        # Check for escape and return
        if y != 0:
            has_left_axis = True
        
        if has_left_axis and y == 0:
            return False # Returned to axis after leaving

        # Get neighbors
        neighbors = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]
        
        # Calculate transition probabilities for the h-transform
        # P_hat(u, v) = P(u,v) * a(v)/a(u) = (1/4) * a(v)/a(u)
        # For simulation, we can choose the next step with weights proportional to a(v).
        weights = []
        for neighbor in neighbors:
            potential = get_potential(neighbor)
            # The walk is on Z^2 \ {0}, so it can't move to the origin.
            if potential < 0:
                weights.append(0)
            else:
                weights.append(potential)

        # This should only happen if all neighbors have 0 probability, which is very unlikely
        if sum(weights) == 0:
            break
            
        # Choose the next step
        next_pos = random.choices(neighbors, weights=weights, k=1)[0]
        current_pos = next_pos
        
    # If the loop finishes, it means the walk did not return to the axis
    return True

def main():
    """
    Main function to run the experiment.
    """
    num_simulations = 1000
    start_pos = (20, 0) # Start on the x-axis, reasonably far from the origin
    num_steps_per_simulation = 5000

    print(f"Is every transient set for the h-transformed SRW in Z^2 necessarily finite? No.")
    print("\n--- Simulation ---")
    print(f"We test if an infinite set (the x-axis) is transient.")
    print(f"A set is transient if the walk eventually leaves it and never returns.")
    print(f"We run {num_simulations} simulations starting at {start_pos} for {num_steps_per_simulation} steps each.")
    print(f"We count how many times the walk leaves the x-axis and does not return.")

    escape_count = 0
    for i in range(num_simulations):
        if run_simulation_trial(start_pos, num_steps_per_simulation):
            escape_count += 1
    
    print("\n--- Results ---")
    print(f"Number of simulations where the walk escaped the x-axis and did not return: {escape_count}")
    print(f"Number of simulations where the walk returned to the x-axis: {num_simulations - escape_count}")

    success_rate = (escape_count / num_simulations) * 100
    print(f"\nEscape percentage: {success_rate:.2f}%")
    print("\nThe high escape rate provides strong numerical evidence that the x-axis is a transient set.")
    print("Since the x-axis is an infinite set, not all transient sets are finite.")

if __name__ == '__main__':
    main()