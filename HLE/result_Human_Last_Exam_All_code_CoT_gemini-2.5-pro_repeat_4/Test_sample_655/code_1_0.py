import math
import numpy as np
from itertools import product

def solve_container_packing():
    """
    Analyzes different container shapes to find the one that packs the most energy.
    It uses a greedy algorithm to place balls, prioritizing larger balls and central placement.
    """

    # Helper function for squared distance, avoids sqrt and floating point issues
    def dist_sq(p1, p2):
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def get_greedy_packing(shape_params):
        """
        Calculates the packing for a given container using a greedy algorithm.
        """
        shape = shape_params['shape']
        
        # Radii of balls
        R_A, R_B = 2.0, 1.0
        
        # Bounding box for generating potential points
        if shape == 'box':
            max_dim = max(shape_params['L'], shape_params['W'], shape_params['H']) / 2
        elif shape == 'sphere':
            max_dim = shape_params['r']
        elif shape == 'cylinder':
            max_dim = max(shape_params['r'], shape_params['h'] / 2)

        # Generate grid points
        step = 0.5
        coords = np.arange(-max_dim, max_dim + step, step)
        grid_points = product(coords, repeat=3)

        # Determine valid center points for each ball type
        P_A_set, P_B_set = set(), set()
        
        for p in grid_points:
            px, py, pz = p
            # Check for Type B (r=1)
            # A ball fits if its center is far enough from all boundaries
            fits_B = False
            if shape == 'box':
                L, W, H = shape_params['L'], shape_params['W'], shape_params['H']
                if abs(px) <= L/2 - R_B and abs(py) <= W/2 - R_B and abs(pz) <= H/2 - R_B:
                    fits_B = True
            elif shape == 'sphere':
                r = shape_params['r']
                if math.sqrt(px**2 + py**2 + pz**2) <= r - R_B:
                    fits_B = True
            elif shape == 'cylinder':
                r, h = shape_params['r'], shape_params['h']
                if math.sqrt(px**2 + py**2) <= r - R_B and abs(pz) <= h/2 - R_B:
                    fits_B = True
            
            if fits_B:
                P_B_set.add(p)
                # Check for Type A (r=2)
                fits_A = False
                if shape == 'box':
                    L, W, H = shape_params['L'], shape_params['W'], shape_params['H']
                    if abs(px) <= L/2 - R_A and abs(py) <= W/2 - R_A and abs(pz) <= H/2 - R_A:
                        fits_A = True
                elif shape == 'sphere':
                    r = shape_params['r']
                    if math.sqrt(px**2 + py**2 + pz**2) <= r - R_A:
                        fits_A = True
                elif shape == 'cylinder':
                    r, h = shape_params['r'], shape_params['h']
                    if math.sqrt(px**2 + py**2) <= r - R_A and abs(pz) <= h/2 - R_A:
                        fits_A = True
                
                if fits_A:
                    P_A_set.add(p)

        # Greedy placement, sorted by distance from origin
        key_func = lambda p: (p[0]**2 + p[1]**2 + p[2]**2, p[0], p[1], p[2])
        sorted_A_points = sorted(list(P_A_set), key=key_func)
        sorted_B_points = sorted(list(P_B_set), key=key_func)

        # Phase 1: Place Type A balls
        placed_A = []
        for p_a in sorted_A_points:
            if all(dist_sq(p_a, placed_p) >= (R_A + R_A)**2 for placed_p in placed_A):
                placed_A.append(p_a)

        # Phase 2: Place Type B balls
        placed_B = []
        for p_b in sorted_B_points:
            # Check conflict with placed Type A balls
            if not all(dist_sq(p_b, p_a) >= (R_A + R_B)**2 for p_a in placed_A):
                continue
            # Check conflict with placed Type B balls
            if all(dist_sq(p_b, placed_p) >= (R_B + R_B)**2 for placed_p in placed_B):
                placed_B.append(p_b)
        
        num_A = len(placed_A)
        num_B = len(placed_B)
        total_energy = 10 * num_A + 1 * num_B
        
        return total_energy, num_A, num_B

    # --- Main Execution ---
    candidates = [
        {'shape': 'sphere', 'r': 9.0},
        {'shape': 'box', 'L': 13.0, 'W': 13.0, 'H': 13.0},
        {'shape': 'cylinder', 'r': 7.5, 'h': 14.5},
        {'shape': 'box', 'L': 10.0, 'W': 10.0, 'H': 21.0}
    ]

    best_energy = -1
    best_config = None

    for params in candidates:
        energy, nA, nB = get_greedy_packing(params)
        
        if energy > best_energy:
            best_energy = energy
            best_config = {
                'params': params,
                'nA': nA,
                'nB': nB,
                'energy': energy
            }

    # Format output string
    p = best_config['params']
    if p['shape'] == 'box':
        C = f"box {p['L']}x{p['W']}x{p['H']}"
    elif p['shape'] == 'sphere':
        C = f"sphere r={p['r']}"
    elif p['shape'] == 'cylinder':
        C = f"cylinder r={p['r']}, h={p['h']}"
    
    a = best_config['nA']
    b = best_config['nB']
    
    result_string = f"[{C}]{a};{b}"
    print(result_string)
    # The final answer format for the platform
    print(f"<<<[{C}]{a};{b}>>>")

solve_container_packing()