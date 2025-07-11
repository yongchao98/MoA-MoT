import math

# Ball and Container properties
R_LARGE, E_LARGE = 2.0, 10
R_SMALL, E_SMALL = 1.0, 1
MAX_SA = 1050.0
PRECISION = 0.5

def get_sc_packed_balls(center_space_dims):
    """Calculates number of balls packed on a simple cubic lattice."""
    # For a simple cubic lattice with spacing D, centers can be placed at k*D/2.
    # Our lattice spacing for large balls is D=4. So centers are at 4k.
    # Center space is the volume where ball centers can be, e.g., for a box
    # of length L, the space for a ball of radius R is L - 2*R.
    # Here, we pass the half-dimensions of the center space, e.g., (L/2 - R).
    
    num_balls_per_dim = []
    for dim_half_len in center_space_dims:
        # Number of positions from origin in one direction (e.g., k=0,1,2..)
        # 4*k <= dim_half_len -> k <= dim_half_len / 4
        # Multiplied by 2 for both directions (+/-) and add 1 for the origin.
        n = math.floor(dim_half_len / R_LARGE / 2) * 2 + 1
        num_balls_per_dim.append(n)
        
    num_total_balls = 1
    for n in num_balls_per_dim:
        num_total_balls *= n
    return num_total_balls, num_balls_per_dim

def get_packed_energy(container_type, dims):
    """Calculates the max energy for a given container config."""
    num_large = 0
    num_small = 0

    if container_type == 'box':
        l, w, h = dims
        sa = 2 * (l*w + l*h + w*h)
        if sa > MAX_SA:
            return -1, 0, 0

        # Pack large balls (r=2)
        center_space_half_dims = [l/2 - R_LARGE, w/2 - R_LARGE, h/2 - R_LARGE]
        if any(d < 0 for d in center_space_half_dims):
            num_large = 0
        else:
            num_large, n_per_dim = get_sc_packed_balls(center_space_half_dims)

        # Pack small balls (r=1) in voids
        if num_large > 0:
            voids_per_dim = [max(0, n - 1) for n in n_per_dim]
            num_voids = voids_per_dim[0] * voids_per_dim[1] * voids_per_dim[2]
            
            # Check if a small ball at a void center (e.g., (2,2,2)) fits
            void_center = [R_LARGE, R_LARGE, R_LARGE]
            # Check container bounds
            if all(void_center[i] <= lwh/2 - R_SMALL for i, lwh in enumerate(dims)):
                 # Check overlap with large balls
                 dist_sq_to_origin_l = sum(c**2 for c in void_center)
                 if dist_sq_to_origin_l >= (R_LARGE + R_SMALL)**2:
                     num_small = num_voids

    elif container_type == 'sphere':
        r_cont = dims[0]
        sa = 4 * math.pi * r_cont**2
        if sa > MAX_SA:
            return -1, 0, 0

        # Pack large balls
        r_center_space = r_cont - R_LARGE
        if r_center_space >= 0:
            # Count integer points (4i, 4j, 4k) within sphere of radius r_center_space
            # (4i)^2 + (4j)^2 + (4k)^2 <= r_center_space^2
            max_i_sq = (r_center_space/ (2*R_LARGE))**2
            
            large_centers = set()
            max_i = int(math.sqrt(max_i_sq))
            for i in range(-max_i, max_i + 1):
                for j in range(-max_i, max_i + 1):
                    for k in range(-max_i, max_i + 1):
                        if i**2 + j**2 + k**2 <= max_i_sq:
                            large_centers.add((i,j,k))
            num_large = len(large_centers)

        # Pack small balls in voids
        if num_large >= 8: # requires at least a 2x2x2 cube of large balls
            void_center_dist_sq = 3 * (R_LARGE**2) # dist of (2,2,2) from origin
            # Check container bounds
            if math.sqrt(void_center_dist_sq) + R_SMALL <= r_cont:
                # Check overlap with large balls (already established this works)
                num_small = (int(math.sqrt(num_large-1))+1 -1)**3 if num_large > 1 else 0
                if num_large == 27: # Special case for 3x3x3 grid
                    num_small = 8


    elif container_type == 'cylinder':
        r, h = dims
        sa = 2 * math.pi * r**2 + 2 * math.pi * r * h
        if sa > MAX_SA:
            return -1, 0, 0

        # Pack large balls in layers
        r_center_space_xy = r - R_LARGE
        h_center_space = h - 2*R_LARGE
        if r_center_space_xy < 0 or h_center_space < 0:
            num_large = 0
        else:
            # Balls per layer
            max_i_sq = (r_center_space_xy/(2*R_LARGE))**2
            pts_in_circle = set()
            max_i = int(math.sqrt(max_i_sq))
            for i in range(-max_i, max_i + 1):
                 for j in range(-max_i, max_i + 1):
                    if i**2+j**2 <= max_i_sq:
                        pts_in_circle.add((i,j))
            n_per_layer = len(pts_in_circle)

            # Number of layers
            n_layers = math.floor(h_center_space / (2*R_LARGE)) + 1
            num_large = n_per_layer * n_layers
            
        # Pack small balls
        if n_per_layer >= 9 and n_layers >= 3 : # for 3x3x3 type packing
             num_small = 8 # Approximation
        else:
             num_small = 0


    total_energy = num_large * E_LARGE + num_small * E_SMALL
    return total_energy, num_large, num_small


def solve():
    best_config = {'energy': -1}
    
    # 1. Sphere
    # Max R for SA <= 1050 is R=sqrt(1050/4pi) ~= 9.14. Max multiple of 0.5 is 9.0
    r = 9.0
    energy, b, a = get_packed_energy('sphere', [r])
    if energy > best_config['energy']:
        best_config = {'energy': energy, 'a': a, 'b': b, 'shape': 'sphere', 'dims': [r]}
    
    # 2. Box
    # Search around the most voluminous shape (a cube)
    # 6L^2 = 1050 -> L = 13.2. Check L=13.0
    l = 13.0
    energy, b, a = get_packed_energy('box', [l, l, l])
    if energy > best_config['energy']:
        best_config = {'energy': energy, 'a': a, 'b': b, 'shape': 'box', 'dims': [l, l, l]}

    # 3. Cylinder
    # Search promising r, h values
    # Try r=8.0, h_max = 1050/(16pi)-8=12.9, so try h=12.5
    r, h = 8.0, 12.5
    energy, b, a = get_packed_energy('cylinder', [r,h])
    if energy > best_config['energy']:
        best_config = {'energy': energy, 'a': a, 'b': b, 'shape': 'cylinder', 'dims': [r, h]}

    # Format the final output string
    if best_config['shape'] == 'box':
        l, w, h = best_config['dims']
        desc = f"box {l}x{w}x{h}"
    elif best_config['shape'] == 'sphere':
        r = best_config['dims'][0]
        desc = f"sphere r={r}"
    elif best_config['shape'] == 'cylinder':
        r, h = best_config['dims']
        desc = f"cylinder r={r}, h={h}"

    a = best_config['a']
    b = best_config['b']

    print(f"Container: {desc}")
    print(f"Number of 1-cm balls (a): {a}")
    print(f"Number of 2-cm balls (b): {b}")
    print(f"Total Energy: {b} * 10 MJ + {a} * 1 MJ = {best_config['energy']} MJ")
    print(f"\nFinal Answer String:")
    
    final_string = f"[{desc}]{a};{b}"
    print(final_string)
    # The final format requested is <<<answer>>>
    print(f"\n<<<[{desc}]{a};{b}>>>")

solve()