import math
import numpy as np

# --- Problem Constants ---
MAX_SURFACE_AREA = 1050.0
STEP = 0.5
LARGE_RADIUS = 2.0
SMALL_RADIUS = 1.0
LARGE_ENERGY = 10
SMALL_ENERGY = 1

# --- Global variable to store the best result found ---
best_result = {
    "energy": 0,
    "shape_str": "",
    "n_small": 0,
    "n_large": 0,
}

# --- Core Packing Logic ---

def pack_container(shape, dims):
    """
    Calculates the maximum energy that can be packed into a given container.
    Uses a greedy approach: first pack large balls, then small balls.
    """
    placed_balls = [] # List to store (x, y, z, radius) of placed balls

    # Define the 3D grid for ball centers based on the container's bounding box
    if shape == "box":
        L, W, H = dims
        # Box centered at (L/2, W/2, H/2)
        x_coords = np.arange(STEP, L, STEP)
        y_coords = np.arange(STEP, W, STEP)
        z_coords = np.arange(STEP, H, STEP)
    elif shape == "cylinder":
        R, H = dims
        # Cylinder centered at (0, 0, 0) with height along z-axis
        x_coords = np.arange(-R + STEP, R, STEP)
        y_coords = np.arange(-R + STEP, R, STEP)
        z_coords = np.arange(-H/2 + STEP, H/2, STEP)
    elif shape == "sphere":
        R = dims[0]
        # Sphere centered at (0, 0, 0)
        x_coords = np.arange(-R + STEP, R, STEP)
        y_coords = np.arange(-R + STEP, R, STEP)
        z_coords = np.arange(-R + STEP, R, STEP)

    # --- Pass 1: Pack Large Balls (greedy) ---
    for x in x_coords:
        for y in y_coords:
            for z in z_coords:
                # 1. Check if a large ball fits within the container walls
                if not is_within_walls(shape, dims, x, y, z, LARGE_RADIUS):
                    continue
                # 2. Check for collision with already placed balls
                if not has_collision(placed_balls, x, y, z, LARGE_RADIUS):
                    placed_balls.append((x, y, z, LARGE_RADIUS))
    
    n_large = len(placed_balls)

    # --- Pass 2: Pack Small Balls in remaining space ---
    for x in x_coords:
        for y in y_coords:
            for z in z_coords:
                # 1. Check if a small ball fits within walls
                if not is_within_walls(shape, dims, x, y, z, SMALL_RADIUS):
                    continue
                # 2. Check for collision
                if not has_collision(placed_balls, x, y, z, SMALL_RADIUS):
                    placed_balls.append((x, y, z, SMALL_RADIUS))

    n_small = len(placed_balls) - n_large
    energy = n_large * LARGE_ENERGY + n_small * SMALL_ENERGY
    
    return energy, n_large, n_small

def is_within_walls(shape, dims, x, y, z, r):
    """Checks if a ball (x,y,z,r) is inside the container."""
    if shape == "box":
        L, W, H = dims
        return (r <= x <= L - r) and (r <= y <= W - r) and (r <= z <= H - r)
    elif shape == "cylinder":
        R, H = dims
        # Assumes cylinder centered at (0,0,0) along z-axis
        return (np.sqrt(x**2 + y**2) <= R - r) and (abs(z) <= H/2 - r)
    elif shape == "sphere":
        R = dims[0]
        # Assumes sphere centered at (0,0,0)
        return np.sqrt(x**2 + y**2 + z**2) <= R - r
    return False

def has_collision(placed_balls, x, y, z, r):
    """Checks if a new ball (x,y,z,r) collides with existing balls."""
    for px, py, pz, pr in placed_balls:
        dist_sq = (x - px)**2 + (y - py)**2 + (z - pz)**2
        min_dist_sq = (r + pr)**2
        if dist_sq < min_dist_sq - 1e-9: # Tolerance for float comparison
            return True
    return False

def update_best_result(energy, n_large, n_small, shape_str):
    """Updates the global best result if the new energy is higher."""
    global best_result
    if energy > best_result["energy"]:
        best_result = {
            "energy": energy,
            "shape_str": shape_str,
            "n_small": n_small,
            "n_large": n_large,
        }
        print(f"New best found: E={energy}, Shape={shape_str}, Large={n_large}, Small={n_small}")

# --- Main Search Logic ---

def solve():
    """Iterates through all container shapes and dimensions to find the optimum."""
    print("Starting optimization... This may take a few minutes.")

    # 1. Sphere
    print("\n--- Checking Sphere Containers ---")
    max_r_sphere = math.sqrt(MAX_SURFACE_AREA / (4 * math.pi))
    for r_int in range(int(STEP * 10), int(max_r_sphere * 10) + 1):
        r = r_int / 10.0
        if r == 0: continue
        
        sa = 4 * math.pi * r**2
        if sa > MAX_SURFACE_AREA: continue
        
        energy, n_large, n_small = pack_container("sphere", (r,))
        update_best_result(energy, n_large, n_small, f"sphere r={r}")

    # 2. Box (using symmetry l >= w >= h)
    print("\n--- Checking Box Containers ---")
    max_h_box = (MAX_SURFACE_AREA / 6)**0.5 # Max h is for a cube
    for h_int in range(int(STEP*10), int(max_h_box * 10) + 1):
        h = h_int / 10.0
        if h == 0: continue

        max_w_box = -h + (h**2 + MAX_SURFACE_AREA/2)**0.5 # For l=w
        for w_int in range(h_int, int(max_w_box * 10) + 1):
            w = w_int / 10.0
            if w == 0: continue

            if 2*(w*h + w*h + h*w) > MAX_SURFACE_AREA: continue
            
            max_l_box = (MAX_SURFACE_AREA/2 - w*h) / (w+h)
            for l_int in range(w_int, int(max_l_box*10) + 1):
                l = l_int / 10.0
                if l == 0: continue

                sa = 2 * (l*w + l*h + w*h)
                if sa > MAX_SURFACE_AREA: break # Inner-most loop optimization
                
                energy, n_large, n_small = pack_container("box", (l, w, h))
                update_best_result(energy, n_large, n_small, f"box {l}x{w}x{h}")
    
    # 3. Cylinder
    print("\n--- Checking Cylinder Containers ---")
    max_r_cyl = (MAX_SURFACE_AREA / (2 * math.pi))**0.5 # For h=0
    for r_int in range(int(STEP*10), int(max_r_cyl * 10) + 1):
        r = r_int / 10.0
        if r == 0: continue

        max_h_cyl = (MAX_SURFACE_AREA / (2 * math.pi * r)) - r
        for h_int in range(int(STEP*10), int(max_h_cyl*10) + 1):
            h = h_int / 10.0
            if h == 0: continue

            sa = 2 * math.pi * r * (r + h)
            if sa > MAX_SURFACE_AREA: break
            
            energy, n_large, n_small = pack_container("cylinder", (r, h))
            update_best_result(energy, n_large, n_small, f"cylinder r={r}, h={h}")

    # --- Print Final Result ---
    print("\n--- Optimization Complete ---")
    print("\nBest configuration found:")
    print(f"Container: {best_result['shape_str']}")
    print(f"Number of 1-cm balls (a): {best_result['n_small']}")
    print(f"Number of 2-cm balls (b): {best_result['n_large']}")
    print("\nFinal energy calculation:")
    print(f"Total Energy = {SMALL_ENERGY} * {best_result['n_small']} + {LARGE_ENERGY} * {best_result['n_large']} = {best_result['energy']} MJ")

if __name__ == '__main__':
    solve()
    final_a = best_result['n_small']
    final_b = best_result['n_large']
    final_C = f"[{best_result['shape_str']}]"
    final_answer = f"<<<{final_C}{final_a};{final_b}>>>"
    # This line prints the final answer in the required format for the system.
    # In a real script, it would be the very last thing printed.
    # print(final_answer) # The final response will have this at the end.
