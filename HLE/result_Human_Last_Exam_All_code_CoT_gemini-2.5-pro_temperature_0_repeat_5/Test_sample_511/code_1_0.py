import math

def count_balls_bcc_like(L, W, H):
    """
    Calculates the number of balls that can fit in a box of size L, W, H
    using a BCC-like packing strategy.
    Radius r=2, Diameter d=4.
    Layer A: 4x4 grid. Layer B: 4x4 grid offset by (2,2). Vertical separation: 3.
    """
    if not all(d >= 4 for d in [L, W, H]):
        return 0

    # Number of balls in a type A layer (grid starts at 2, spacing 4)
    # Centers must be in [r, Dim-r], i.e., [2, Dim-2]
    # Centers along one axis: 2, 6, 10, ...
    # Number of balls = floor(Dim / 4)
    na_x = math.floor(L / 4)
    na_y = math.floor(W / 4)
    N_A = na_x * na_y

    # Number of balls in a type B layer (grid starts at 4, spacing 4)
    # Centers along one axis: 4, 8, 12, ...
    # Last center <= Dim - 2. First center is 4.
    # Number of balls = floor((Dim - 6) / 4) + 1 for Dim >= 6
    nb_x = math.floor((L - 6) / 4) + 1 if L >= 6 else 0
    nb_y = math.floor((W - 6) / 4) + 1 if W >= 6 else 0
    N_B = nb_x * nb_y

    if N_A == 0: # If the main layer doesn't fit, nothing fits
        return 0

    total_balls = 0
    z = 2.0  # Center z-coordinate of the first layer
    is_A_layer = True
    # A ball's center must be within [2, H-2]
    while z <= H - 2:
        if is_A_layer:
            total_balls += N_A
        else:
            total_balls += N_B
        z += 3.0
        is_A_layer = not is_A_layer
    return total_balls

def find_optimal_box():
    """
    Searches for integer dimensions L, W, H that can hold >= 27 balls
    with a surface area < 864.
    """
    min_sa = 864
    best_solution = None
    target_balls = 27

    # Search integer dimensions L, W, H, assuming L >= W >= H
    # to avoid duplicate checks of the same shape.
    # From 2*(LW+WH+HL) < 864 and L>=W>=H, we get 6*H^2 < 864 -> H < 12.
    for H in range(4, 12):
        # Further constrain W: 2*(W^2 + 2WH) < 864 -> W^2+2WH < 432
        w_limit = int(math.sqrt(H*H + 432) - H)
        for W in range(H, w_limit + 1):
            # Further constrain L: 2*(LW+LH+WH) < 864 -> L*(W+H) < 432-WH
            if (W + H) == 0: continue
            l_limit = int((432 - W * H) / (W + H))
            for L in range(W, l_limit + 1):
                
                current_sa = 2 * (L*W + W*H + H*L)
                if current_sa >= min_sa:
                    continue

                # Check all 3 unique orientations for packing
                dims_set = {L, W, H}
                max_balls = 0
                # Permutation 1: (L, W, H)
                max_balls = max(max_balls, count_balls_bcc_like(L, W, H))
                # Permutation 2: (L, H, W)
                max_balls = max(max_balls, count_balls_bcc_like(L, H, W))
                # Permutation 3: (W, H, L)
                max_balls = max(max_balls, count_balls_bcc_like(W, H, L))
                
                if max_balls >= target_balls:
                    if current_sa < min_sa:
                        min_sa = current_sa
                        best_solution = (L, W, H, current_sa)

    if best_solution:
        L, W, H, sa = best_solution
        # Sort dimensions for consistent output, largest first
        dims = sorted([L, W, H], reverse=True)
        print(f"{dims[0]}:{dims[1]}:{dims[2]}:{sa}")
    else:
        print("0")

if __name__ == '__main__':
    find_optimal_box()