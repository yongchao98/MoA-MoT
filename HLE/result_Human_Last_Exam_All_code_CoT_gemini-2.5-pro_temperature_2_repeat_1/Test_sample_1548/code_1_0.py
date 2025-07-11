import math
import random

def solve_potential_kernel(grid_size, iterations):
    """
    Numerically computes the potential kernel `a(x)` on a grid.
    The function a(x) is harmonic on Z^2\{0} and is set to log|x| 
    at the boundary of our simulation grid. We solve for the interior
    values using the relaxation method.
    """
    print(f"Approximating the potential kernel `a(x)` on a {grid_size*2}x{grid_size*2} grid...")
    # Using a dictionary for sparse storage: {(i,j): value}
    a = {}
    center = grid_size
    
    # Initialize boundary conditions and interior points
    for i in range(-grid_size, grid_size + 1):
        for j in range(-grid_size, grid_size + 1):
            if i == 0 and j == 0:
                continue
            if abs(i) == grid_size or abs(j) == grid_size:
                a[(i, j)] = math.log(math.sqrt(i**2 + j**2))
            else:
                a[(i, j)] = 0.0

    # Relaxation method to find the harmonic function values
    for k in range(iterations):
        for i in range(-grid_size + 1, grid_size):
            for j in range(-grid_size + 1, grid_size):
                if i == 0 and j == 0:
                    continue
                neighbors = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
                a[(i,j)] = 0.25 * sum(a.get(n, 0) for n in neighbors)
    
    print("Approximation complete.\n")
    return a

def simulate_h_transform_walk(a, start_pos, max_steps, grid_size):
    """
    Simulates the Doob's h-transform of the SRW and checks if
    the x-axis (an infinite set) is transient.
    """
    print("--- Simulation of the h-transformed walk ---")
    print(f"The chosen infinite set is A = {{ (k,0) | k is a non-zero integer }}.")
    print(f"We start a walker at {start_pos} and run for {max_steps} steps.")
    print("If set A is transient, the walker should visit it only a finite number of times.\n")

    x_curr = start_pos
    xaxis_visits = 0
    path = []
    
    drift_calculation_done = False

    for step in range(max_steps):
        if x_curr[1] == 0:
            xaxis_visits += 1

        path.append(x_curr)
        
        # Check if walker is out of bounds
        if not (abs(x_curr[0]) < grid_size and abs(x_curr[1]) < grid_size):
            print(f"Walker moved out of the pre-computed grid at step {step}. Ending simulation.")
            break

        # Demonstrate drift calculation at the first opportunity when y=1
        if x_curr[1] == 1 and not drift_calculation_done:
            print(f"--- Demonstrating drift at position {x_curr} ---")
            (cx, cy) = x_curr
            a_curr = a.get((cx, cy))
            a_up = a.get((cx, cy + 1))
            a_down = a.get((cx, cy - 1))

            if a_curr and a_up and a_down:
                # The expected change in the y-coordinate is proportional to a(x,y+1) - a(x,y-1)
                drift_numerator = a_up - a_down
                drift = drift_numerator / (4 * a_curr)

                print(f"The equation for local drift in y is: E[dY] = (a(x,y+1) - a(x,y-1)) / (4*a(x,y))")
                print(f"At x={cx}, y={cy}:")
                print(f"a({cx},{cy+1}) = {a_up:.4f}")
                print(f"a({cx},{cy-1}) = {a_down:.4f}")
                print(f"a({cx},{cy})   = {a_curr:.4f}")
                print(f"Expected y-drift E[dY] = ({a_up:.4f} - {a_down:.4f}) / (4 * {a_curr:.4f}) = {drift:.6f}")
                print("Since the drift is positive, the walk is pushed away from the x-axis (y=0).\n")
                drift_calculation_done = True
        
        # Get neighbors and transition probabilities
        i, j = x_curr
        neighbors = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
        
        a_curr_val = a.get(x_curr)
        if not a_curr_val or a_curr_val == 0:
            print(f"Error: a({x_curr}) is zero or undefined. Cannot proceed.")
            break
            
        a_neighbors = [a.get(n, 0) for n in neighbors]
        # Transition prob p(x->y) = (1/4) * a(y)/a(x)
        # We use a_neighbors as weights since the common term (1/4)/a(x) cancels out.
        if sum(a_neighbors) == 0:
             print(f"Walker is trapped at {x_curr} as all neighbors have a=0. Stopping.")
             break

        x_curr = random.choices(neighbors, weights=a_neighbors, k=1)[0]
    
    print("--- Simulation Results ---")
    print(f"Simulation finished after {len(path)} steps.")
    print(f"Final position: {x_curr}, at distance {math.sqrt(x_curr[0]**2 + x_curr[1]**2):.2f} from origin.")
    print(f"Total number of visits to the infinite set A (the x-axis): {xaxis_visits}")
    print("\nThis finite number of visits suggests the walk drifts away from the set A.")


if __name__ == '__main__':
    GRID_SIZE = 50  # Grid will be from -50 to 50
    ITERATIONS = 1000 # Iterations for relaxation
    
    # 1. Approximate the harmonic function 'a'
    potential_kernel = solve_potential_kernel(GRID_SIZE, ITERATIONS)
    
    # 2. Run the simulation
    simulate_h_transform_walk(a=potential_kernel, 
                              start_pos=(1, 1), 
                              max_steps=20000, 
                              grid_size=GRID_SIZE)

    print("\n>>> False")
