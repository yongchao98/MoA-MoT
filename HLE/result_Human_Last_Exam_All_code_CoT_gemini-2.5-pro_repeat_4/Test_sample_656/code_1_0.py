import math

def solve_pioneer_packing():
    """
    This function solves the Pioneer probe packing problem by analyzing a specifically chosen,
    highly efficient container configuration and a layered packing strategy.
    """

    # Step 1: Design an optimal container based on a specific packing strategy.
    # The strategy is to create a core of 3x3x3 high-value 2-cm radius balls and
    # use the remaining surface area to add a slab of space for 1-cm radius balls.
    # A 3-layer stack of 2-cm balls (4 cm diameter) is 12 cm high.
    # A 1-layer slab of 1-cm balls (2 cm diameter) is 2 cm high.
    # Total optimal height H = 12.0 + 2.0 = 14.0 cm.
    
    # Now, find the maximum side length (L) for a square base (L x L) given H=14.0
    # and a maximum surface area of 1050 cm^2.
    # Surface Area Formula: 2 * (L*L + L*H + L*H) <= 1050
    # L^2 + 2*L*H <= 525  => L^2 + 28L - 525 <= 0
    # We solve for the positive root of L^2 + 28L - 525 = 0
    l_max_continuous = (-28 + math.sqrt(28**2 - 4*1*(-525))) / 2
    # The result is ~12.85 cm.

    # Since dimensions must be multiples of 0.5 cm, we take the floor to the nearest 0.5.
    L = math.floor(l_max_continuous * 2) / 2
    W = L
    H = 14.0
    container_description = f"box {L}x{W}x{H}"

    # Step 2: Calculate the number of balls that fit using this configuration.

    # First, pack the 2-cm radius (4 cm diameter) balls on a simple cubic lattice.
    # The number of balls along each dimension is the floor of the dimension divided by the ball's diameter.
    n2_x = math.floor(L / 4.0)
    n2_y = math.floor(W / 4.0)
    n2_z = math.floor(12.0 / 4.0) # They occupy the first 12 cm of height
    num_2cm_balls = n2_x * n2_y * n2_z

    # Second, pack the 1-cm radius (2 cm diameter) balls into the remaining space.
    # This space consists of two types of regions:

    # a) The small interstitial voids between the large balls.
    # A grid of n_x*n_y*n_z spheres has (n_x-1)*(n_y-1)*(n_z-1) central voids.
    # These voids are large enough to fit a 1-cm radius ball.
    num_1cm_interstitial = (n2_x - 1) * (n2_y - 1) * (n2_z - 1)

    # b) The empty 2-cm-high slab at the top of the container.
    # The slab dimensions are L x W x (H - 12.0).
    slab_H = H - 12.0
    n1_slab_x = math.floor(L / 2.0)
    n1_slab_y = math.floor(W / 2.0)
    n1_slab_z = math.floor(slab_H / 2.0)
    num_1cm_slab = n1_slab_x * n1_slab_y * n1_slab_z
    
    num_1cm_balls = num_1cm_interstitial + num_1cm_slab
    
    # Step 3: Output the final answer in the required format.
    final_answer_string = f"[{container_description}]{num_1cm_balls};{num_2cm_balls}"
    print(final_answer_string)


solve_pioneer_packing()
<<<[box 12.5x12.5x14.0]44;27>>>