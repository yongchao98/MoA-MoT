import math

def solve_pioneer_packing():
    """
    This function analyzes several candidate container designs to find the one
    that maximizes the total energy packed, subject to a surface area constraint of 1050 cm^2.
    The analysis prioritizes packing the high-energy 2-cm radius balls.
    """

    # --- Properties ---
    # Large ball: radius = 2 cm, diameter = 4 cm, energy = 20 MJ
    # Small ball: radius = 1 cm, diameter = 2 cm, energy = 1 MJ
    # Surface Area <= 1050 cm^2
    # Dimensions are multiples of 0.5 cm

    best_config = {
        "description": "",
        "num_small": 0,
        "num_large": 0,
        "energy": 0
    }

    # --- Candidate 1: Box optimized for grid packing of large balls ---
    # A box with dimensions that are multiples of the large ball's diameter (4 cm).
    # Let's test a 16x16x8 cm box.
    l, w, h = 16.0, 16.0, 8.0
    sa1 = 2 * (l*w + l*h + w*h)
    if sa1 <= 1050:
        # Number of large balls is determined by grid packing.
        nl1 = (l / 4) * (w / 4) * (h / 4)
        # There is no leftover space for small balls.
        ns1 = 0
        energy1 = nl1 * 20 + ns1 * 1
        if energy1 > best_config["energy"]:
            best_config.update({
                "description": f"box {l}x{w}x{h}",
                "num_small": int(ns1),
                "num_large": int(nl1),
                "energy": energy1
            })

    # --- Candidate 2: Cylinder optimized for dense hexagonal packing ---
    # A known efficient packing fits 7 spheres (radius r) in a cylinder of radius 3r.
    # For large balls (r=2 cm), this requires a container cylinder of radius R=6 cm.
    # We can stack layers, where each layer has a height of one ball diameter (4 cm).
    # We find the maximum number of layers (N) allowed by the surface area constraint.
    # SA = 2*pi*R^2 + 2*pi*R*H, where H = N * 4.
    # 1050 >= 2*pi*6^2 + 2*pi*6*(4*N)
    # 1050 >= 72*pi + 48*pi*N
    # (1050 - 226.2) / 150.8 >= N  -->  5.46 >= N
    # So we can have a maximum of 5 layers.
    cyl_r, cyl_h = 6.0, 5 * 4.0 # 5 layers * 4cm/layer = 20cm height
    sa2 = 2 * math.pi * cyl_r**2 + 2 * math.pi * cyl_r * cyl_h
    if sa2 <= 1050:
        # Number of large balls = 7 balls/layer * 5 layers.
        nl2 = 7 * 5
        # This packing is very tight, leaving no significant space for small balls.
        ns2 = 0
        energy2 = nl2 * 20 + ns2 * 1
        if energy2 > best_config["energy"]:
            best_config.update({
                "description": f"cylinder r={cyl_r}, h={cyl_h}",
                "num_small": int(ns2),
                "num_large": int(nl2),
                "energy": energy2
            })
    
    # --- Candidate 3: Box for mixed packing (large balls + layer of small balls) ---
    # Let's start with a 12x12x12 cm space to fit a 3x3x3 grid of 27 large balls.
    # This base has SA = 6*12^2 = 864 cm^2. We can extend one dimension.
    # Let's try a 12x12x14 cm box.
    l, w, h = 12.0, 12.0, 14.0
    sa3 = 2 * (l*w + l*h + w*h)
    if sa3 <= 1050:
        # It fits 27 large balls in a 12x12x12 section.
        nl3 = 27
        # This leaves a 12x12x2 cm slab for small balls (diameter 2 cm).
        ns3 = (12 / 2) * (12 / 2) * (2 / 2)
        energy3 = nl3 * 20 + ns3 * 1
        if energy3 > best_config["energy"]:
            best_config.update({
                "description": f"box {l}x{w}x{h}",
                "num_small": int(ns3),
                "num_large": int(nl3),
                "energy": energy3
            })

    # Print the best result found.
    # The format is [Description]a;b where a is num_small_balls, b is num_large_balls
    final_desc = best_config["description"]
    final_ns = best_config["num_small"]
    final_nl = best_config["num_large"]

    print(f"[{final_desc}]{final_ns};{final_nl}")

solve_pioneer_packing()