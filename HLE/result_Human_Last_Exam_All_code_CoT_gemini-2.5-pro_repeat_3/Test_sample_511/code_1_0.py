import math

def solve_and_print_answer():
    """
    This function searches for a more efficient container for the energy balls.
    It establishes a baseline with the initial 12x12x12 box and then searches
    for alternative integer dimensions that offer a smaller surface area while holding
    at least the same number of spheres.
    """
    # --- Problem Parameters ---
    initial_dims = (12, 12, 12)
    sphere_radius = 2.0
    diameter = 2 * sphere_radius

    # --- Baseline Calculation ---
    l_init, w_init, h_init = initial_dims
    sa_initial = 2 * (l_init*w_init + w_init*h_init + h_init*l_init)
    
    # Calculate the number of spheres the initial box can hold with simple cubic packing
    n_l_init = math.floor((l_init - diameter) / diameter) + 1 if l_init >= diameter else 0
    n_w_init = math.floor((w_init - diameter) / diameter) + 1 if w_init >= diameter else 0
    n_h_init = math.floor((h_init - diameter) / diameter) + 1 if h_init >= diameter else 0
    num_spheres_required = n_l_init * n_w_init * n_h_init

    # --- Optimization Search ---
    best_box = None
    min_sa = sa_initial

    # The search space for integer dimensions is centered around the initial values.
    # We are looking for a box with SA < sa_initial, which implies its volume
    # will also be smaller than the initial volume (1728).
    # Smallest dimension must be at least the sphere diameter (4cm).
    # We can set a reasonable upper bound for the search.
    search_limit = 20 

    for l in range(4, search_limit + 1):
      # To avoid duplicate permutations like (10,11,12) and (11,10,12),
      # we enforce l <= w <= h
      for w in range(l, search_limit + 1):
        for h in range(w, search_limit + 1):
          sa_current = 2 * (l*w + w*h + h*l)
          
          # We only care about boxes with a smaller surface area
          if sa_current < min_sa:
            # Check if this box can hold enough spheres with simple cubic packing
            nl_cand = math.floor((l - diameter) / diameter) + 1 if l >= diameter else 0
            nw_cand = math.floor((w - diameter) / diameter) + 1 if w >= diameter else 0
            nh_cand = math.floor((h - diameter) / diameter) + 1 if h >= diameter else 0
            
            num_spheres_cand = nl_cand * nw_cand * nh_cand

            if num_spheres_cand >= num_spheres_required:
              # Found a better box
              min_sa = sa_current
              best_box = (l, w, h)

    # --- Final Output ---
    if best_box is None:
      # No better box was found in the search
      print("0")
    else:
      # A better box was found, format the output
      # The problem asks for a:b:c:d, so we provide the found dimensions.
      dims = sorted(list(best_box), reverse=True)
      l_best, w_best, h_best = dims[0], dims[1], dims[2]
      sa_best = int(round(min_sa))
      print(f"{l_best}:{w_best}:{h_best}:{sa_best}")

solve_and_print_answer()
<<<0>>>