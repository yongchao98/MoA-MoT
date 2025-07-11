import math

def solve_pioneer_probe_packing():
    """
    Solves the energy packing problem by searching for the optimal container
    dimensions for both a box and a cylinder, maximizing the total energy.
    """

    # The objective is to maximize energy E = 1*n1 + 20*n2.
    # The problem is simplified by considering the bounding box volumes for the spheres:
    # - 1-cm radius ball: 2x2x2 cm cube, volume = 8 cm^3, energy = 1 MJ.
    # - 2-cm radius ball: 4x4x4 cm cube, volume = 64 cm^3, energy = 20 MJ.
    # Energy density of large ball (20/64) > small ball (1/8), so we prioritize large balls.
    
    def calculate_energy_from_volume(volume):
        """
        Calculates the maximum energy that can be packed into a given volume.
        It prioritizes the more energy-dense 2-cm balls.
        """
        if volume < 0:
            return 0, 0, 0
            
        n2 = math.floor(volume / 64)
        remaining_volume = volume - n2 * 64
        n1 = math.floor(remaining_volume / 8)
        energy = n1 * 1 + n2 * 20
        return energy, n1, n2

    max_surface_area = 1050.0
    step = 0.5

    # --- Part 1: Box Container Search ---
    best_box = {'energy': 0}
    
    # Iterate through possible dimensions (multiples of 0.5).
    # To avoid duplicate shapes and reduce search space, we enforce L >= W >= H.
    # A cube gives max volume for SA, so 6*H^2 <= 1050 -> H <= sqrt(175) ~= 13.2.
    # So, we test H up to 13.0 cm.
    h_limit_i = int(13.0 / step)
    for h_i in range(1, h_limit_i + 1):
        H = h_i * step
        
        # Upper bound for W: 4*W*H <= 1050 -> W <= 1050 / (4*H)
        w_limit_i = int((max_surface_area / (4 * H)) / step)
        for w_i in range(h_i, w_limit_i + 1):
            W = w_i * step
            
            # Upper bound for L: L*(2W + 2H) <= 1050 - 2WH
            if (2 * W + 2 * H) == 0: continue
            l_limit = (max_surface_area - 2 * W * H) / (2 * W + 2 * H)
            l_limit_i = int(l_limit / step)
            for l_i in range(w_i, l_limit_i + 1):
                L = l_i * step
                
                surface_area = 2 * (L * W + L * H + W * H)
                if surface_area <= max_surface_area:
                    volume = L * W * H
                    energy, n1, n2 = calculate_energy_from_volume(volume)
                    if energy > best_box['energy']:
                        best_box = {'energy': energy, 'L': L, 'W': W, 'H': H, 'n1': n1, 'n2': n2}

    # --- Part 2: Cylinder Container Search ---
    best_cylinder = {'energy': 0}
    pi = math.pi
    
    # Iterate through possible dimensions for the cylinder.
    # 2*pi*R^2 <= 1050 -> R <= sqrt(1050/(2*pi)) ~= 12.9.
    # We test R up to 12.5 cm.
    r_limit_i = int(12.5 / step)
    for r_i in range(1, r_limit_i + 1):
        R = r_i * step
        
        # Upper bound for H: H <= 1050/(2*pi*R) - R
        if (2 * pi * R) == 0: continue
        h_limit = (max_surface_area / (2 * pi * R)) - R
        if h_limit < step: continue
        h_limit_i = int(h_limit / step)
        
        for h_i in range(1, h_limit_i + 1):
            H = h_i * step
            surface_area = 2 * pi * R**2 + 2 * pi * R * H
            if surface_area <= max_surface_area:
                volume = pi * R**2 * H
                energy, n1, n2 = calculate_energy_from_volume(volume)
                if energy > best_cylinder['energy']:
                    best_cylinder = {'energy': energy, 'R': R, 'H': H, 'n1': n1, 'n2': n2}

    # --- Part 3: Compare and Print the Final Answer ---
    if best_box.get('energy', 0) > best_cylinder.get('energy', 0):
        dims = sorted([best_box['L'], best_box['W'], best_box['H']], reverse=True)
        desc = f"box {dims[0]}x{dims[1]}x{dims[2]}"
        n1, n2 = best_box['n1'], best_box['n2']
        print(f"[{desc}]{n1};{n2}")
    else:
        desc = f"cylinder r={best_cylinder['R']}, h={best_cylinder['H']}"
        n1, n2 = best_cylinder['n1'], best_cylinder['n2']
        print(f"[{desc}]{n1};{n2}")

solve_pioneer_probe_packing()