import math
import random

def simulate_h_transform_walk(steps, initial_pos, test_set_func):
    """
    Simulates Doob's h-transform of a 2D simple random walk.

    This simulation uses an approximation of the potential kernel.
    The true potential kernel a(x) is approximated by h(x) = log(|x|).

    Args:
        steps (int): The number of steps to simulate.
        initial_pos (list): The starting position [x, y].
        test_set_func (function): A function that returns True if a position
                                  is in the set we are testing, False otherwise.
    """
    pos = list(initial_pos)
    visit_count = 0

    # The h-function h(x) for the transform is the potential kernel a(x). For
    # large distances, a(x) behaves like log(|x|). We use this approximation.
    # To handle positions near the origin, we use h(x,y) = log(x^2+y^2).
    # The walk must start at a position with a norm greater than 1.
    def h(p):
        norm_sq = p[0]**2 + p[1]**2
        # The h-function is positive for norm > 1.
        if norm_sq <= 1:
            return 0
        # The potential is proportional to the log of the norm.
        return math.log(norm_sq)

    # Initial check of the starting position
    if h(pos) == 0:
      print("Error: Start position must have a distance > 1 from the origin.")
      return

    print(f"Starting simulation for {steps} steps from position {pos}.")
    print("We are testing if the positive x-axis is a transient set.")
    print("If the set were transient, the number of visits would stay bounded.")
    print("We expect the number of visits to grow, suggesting the set is recurrent.")
    print("-" * 55)
    print("Step         | Distance from Origin | Visits to Set")
    print("-" * 55)

    # Main simulation loop
    for step in range(1, steps + 1):
        if test_set_func(pos):
            visit_count += 1

        # Get current position and neighbors
        x, y = pos
        neighbors = [[x + 1, y], [x - 1, y], [x, y + 1], [x, y - 1]]

        # For Doob's h-transform, the transition probability P(x->y) is proportional to h(y).
        # We calculate the weights for each neighbor based on the h-function.
        h_vals = [h(p) for p in neighbors]
        
        # We create a probability distribution by normalizing the weights.
        # This is a valid Markov chain that is biased towards larger h-values.
        total_h = sum(h_vals)
        if total_h <= 0:
             # This might happen if all neighbors are inside the unit circle.
             # In this case, we move to the neighbor with the largest h-value,
             # effectively pushing the walk away from the origin.
             max_h_val = -1
             next_pos = pos
             for i, h_val in enumerate(h_vals):
                 if h_val > max_h_val:
                     max_h_val = h_val
                     next_pos = neighbors[i]
             pos = next_pos
        else:
            probs = [h_val / total_h for h_val in h_vals]
            # Choose the next position based on the calculated probabilities
            pos = random.choices(neighbors, weights=probs, k=1)[0]
        
        # Periodically print progress
        if step % (steps // 10) == 0:
            distance = math.sqrt(pos[0]**2 + pos[1]**2)
            print(f"{step:<12} | {distance:<20.2f} | {visit_count:<10}")

    print("-" * 55)
    print(f"Simulation finished after {steps} steps.")
    print(f"Final position: [{pos[0]}, {pos[1]}]")
    print(f"Total visits to the positive x-axis: {visit_count}")

def is_on_positive_x_axis(pos):
    """Checks if a position is on the positive x-axis."""
    return pos[1] == 0 and pos[0] > 0

# --- Main execution ---
if __name__ == '__main__':
    # A large number of steps is needed to observe long-term behavior.
    total_steps = 200000
    # Start sufficiently far from the origin so that h(x) ~ log|x| is a good approximation.
    start_pos = [10, 0]

    simulate_h_transform_walk(total_steps, start_pos, is_on_positive_x_axis)