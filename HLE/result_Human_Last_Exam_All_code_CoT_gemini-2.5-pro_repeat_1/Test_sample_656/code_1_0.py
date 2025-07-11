import numpy as np

def solve():
    """
    Solves the energy ball packing problem by searching through possible container
    dimensions and using a lattice-based heuristic to estimate the packed energy.
    """
    
    # The overall plan is to perform an exhaustive search over possible container
    # dimensions, constrained by the maximum surface area. For each valid container,
    # we estimate the maximum energy it can hold using a simplified packing model.
    # This model is a heuristic because solving the true 3D sphere packing problem is
    # computationally very difficult.
    #
    # Step-by-step plan:
    # 1. Iterate through container types: first box, then cylinder.
    # 2. For each type, iterate through all possible dimensions that are multiples of 0.5 cm.
    #    The search range for dimensions is pruned by the surface area constraint (<= 1050 cm^2)
    #    to make the search feasible.
    # 3. For each valid container, a packing heuristic is applied:
    #    a. We assume the balls are placed on a simple cubic lattice. This is a strong
    #       heuristic given the constraint that ball centers must be on a 0.5 cm grid.
    #    b. We prioritize packing the high-energy 2-cm balls first. The heuristic calculates
    #       the largest (nx, ny, nz) grid of 2-cm balls (diameter 4 cm) that can fit.
    #    c. We then pack 1-cm balls in the main interstitial sites of the 2-cm ball lattice.
    #       For a cubic lattice of 2-cm balls, a 1-cm ball can fit at the center of the
    #       4x4x4 cm cube between the larger balls. The number of such sites is (nx-1)*(ny-1)*(nz-1).
    # 4. The total energy (Energy = 20 * n2 + 1 * n1) is calculated for each configuration.
    # 5. The configuration (container type, dimensions, and ball counts) that gives the
    #    highest energy is stored.
    # 6. Finally, the best result found is printed in the required format "[C]a;b".

    max_surface_area = 1050.0
    step = 0.5
    
    best_energy = -1.0
    best_config = {}

    def get_box_packing(l, w, h):
        """Heuristic for packing a box based on a simple cubic lattice."""
        if l < 4.0 or w < 4.0 or h < 4.0:
            nx, ny, nz = 0, 0, 0
        else:
            # A grid of nx balls with diameter 4 needs a length of 4*nx
            nx = int(l / 4.0)
            ny = int(w / 4.0)
            nz = int(h / 4.0)
        
        n2 = nx * ny * nz
        
        # Pack 1-cm balls in the body-centered interstitial sites of the 2-cm ball lattice
        n1 = 0
        if nx > 1 and ny > 1 and nz > 1:
            n1 = (nx - 1) * (ny - 1) * (nz - 1)
            
        energy = 20.0 * n2 + 1.0 * n1
        return energy, n1, n2

    def get_cylinder_packing(r, h):
        """Heuristic for packing a cylinder by fitting a rectangular lattice inside."""
        if r < 2.0 or h < 4.0:
            return 0.0, 0, 0

        nz = int(h / 4.0)
        if nz == 0:
            return 0.0, 0, 0

        max_n2_layer = 0
        nx, ny = 0, 0
        
        # To fit a rectangular grid of 4x4 squares (for 2cm balls), a box of size
        # 4*nx by 4*ny must fit in a circle of radius r.
        # The condition is (2*nx)^2 + (2*ny)^2 <= r^2.
        # We search for the nx, ny pair that maximizes nx*ny.
        max_nx_check = int(r / 2.0)
        for nx_try in range(1, max_nx_check + 1):
            if 4.0 * nx_try**2 > r**2: continue
            max_ny_sq = r**2 / 4.0 - nx_try**2
            if max_ny_sq >= 1.0:
                ny_try = int(np.sqrt(max_ny_sq))
                if nx_try * ny_try > max_n2_layer:
                    max_n2_layer = nx_try * ny_try
                    nx, ny = nx_try, ny_try
        
        n2 = max_n2_layer * nz

        n1 = 0
        if nx > 1 and ny > 1 and nz > 1:
             n1 = (nx - 1) * (ny - 1) * (nz - 1)
        
        energy = 20.0 * n2 + 1.0 * n1
        return energy, n1, n2

    # --- Search for Best Box (L >= W >= H) ---
    # To optimize search, we set L>=W>=H to avoid redundant checks of same box with different dimension orders.
    h_max = int(np.sqrt(max_surface_area / 6.0) / step)
    for h_i in range(int(4.0/step), h_max + 3):
        h = h_i * step
        w_max_bound = np.sqrt(h**2 + max_surface_area) - h
        for w_i in range(h_i, int(w_max_bound / step) + 3):
            w = w_i * step
            if 2 * (w*h) > max_surface_area: continue
            
            l_max_num = max_surface_area / 2.0 - w * h
            l_max_den = w + h
            if l_max_den <= 0: continue
            for l_i in range(w_i, int(l_max_num / l_max_den / step) + 3):
                l = l_i * step
                sa = 2 * (l*w + w*h + h*l)
                if sa > max_surface_area: break 

                energy, n1, n2 = get_box_packing(l, w, h)
                if energy > best_energy:
                    best_energy = energy
                    best_config = {"type": "box", "dims": (l, w, h), "n1": n1, "n2": n2}

    # --- Search for Best Cylinder ---
    r_max = int(np.sqrt(max_surface_area / (2 * np.pi)) / step)
    for r_i in range(int(2.0/step), r_max + 3):
        r = r_i * step
        sa_base = 2 * np.pi * r**2
        if sa_base > max_surface_area: break
        
        h_max_den = 2 * np.pi * r
        if h_max_den <= 0: continue
        h_max_num = max_surface_area - sa_base
        for h_i in range(int(4.0/step), int(h_max_num / h_max_den / step) + 3):
            h = h_i * step
            sa = sa_base + 2 * np.pi * r * h
            if sa > max_surface_area: break
                
            energy, n1, n2 = get_cylinder_packing(r, h)
            if energy > best_energy:
                best_energy = energy
                best_config = {"type": "cylinder", "dims": (r, h), "n1": n1, "n2": n2}
    
    # --- Format and Print Final Answer ---
    if not best_config:
        print("[0]")
        return
        
    cfg = best_config
    n1 = cfg["n1"]
    n2 = cfg["n2"]

    desc = ""
    if cfg["type"] == "box":
        l, w, h = cfg["dims"]
        l_str = str(l) if l != int(l) else str(int(l))
        w_str = str(w) if w != int(w) else str(int(w))
        h_str = str(h) if h != int(h) else str(int(h))
        desc = f"box {l_str}x{w_str}x{h_str}"
    elif cfg["type"] == "cylinder":
        r, h = cfg["dims"]
        r_str = str(r) if r != int(r) else str(int(r))
        h_str = str(h) if h != int(h) else str(int(h))
        desc = f"cylinder r={r_str}, h={h_str}"
    
    final_answer = f"[{desc}]{n1};{n2}"
    print(final_answer)

if __name__ == '__main__':
    solve()
<<<[box 13x13x13]8;27>>>