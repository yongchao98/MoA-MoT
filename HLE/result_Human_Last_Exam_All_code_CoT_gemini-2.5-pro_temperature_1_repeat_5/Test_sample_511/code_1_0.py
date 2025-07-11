import math

def solve():
    """
    Analyzes the sphere packing problem to find a more efficient container.

    The function calculates the properties of the initial container and then
    searches for a new container configuration (defined by the number of balls
    along each axis) that can hold at least the same number of balls but
    with a smaller surface area.

    The packing is assumed to be a simple cubic lattice, as implied by the
    problem's constraints on coordinates and manufacturing precision.
    """
    
    # --- Initial Configuration ---
    L0, W0, H0 = 12, 12, 12  # Initial box dimensions in cm
    R_BALL = 2              # Ball radius in cm
    D_BALL = 2 * R_BALL     # Ball diameter in cm

    # Calculate initial surface area and ball capacity
    SA0 = 2 * (L0*W0 + L0*H0 + W0*H0)
    
    # Number of balls in the initial box using simple cubic packing
    nx0 = math.floor((L0 - D_BALL) / D_BALL) + 1
    ny0 = math.floor((W0 - D_BALL) / D_BALL) + 1
    nz0 = math.floor((H0 - D_BALL) / D_BALL) + 1
    N0 = nx0 * ny0 * nz0
    
    best_config = None
    min_found_sa = float('inf')

    # Search for a better configuration (nx, ny, nz)
    # The search space for N is bounded; as N increases, SA invariably increases.
    # A search limit of N0*2 is more than sufficient.
    for n_balls in range(N0, N0 * 2):
        
        # For each n_balls, find the integer factors (nx, ny, nz) that are
        # closest to each other (most cube-like) to minimize surface area.
        min_sa_term = float('inf')
        best_factors_for_n = None

        # Efficiently find factors
        for nx in range(1, int(n_balls**(1/3.0)) + 2):
            if n_balls % nx == 0:
                rem = n_balls // nx
                for ny in range(nx, int(rem**0.5) + 2):
                    if rem % ny == 0:
                        nz = rem // ny
                        
                        sa_term = nx*ny + nx*nz + ny*nz
                        if sa_term < min_sa_term:
                            min_sa_term = sa_term
                            best_factors_for_n = (nx, ny, nz)
        
        # If a valid factorization was found, calculate the SA
        if best_factors_for_n:
            nx, ny, nz = best_factors_for_n
            
            # The minimum box dimensions required are integers
            L = D_BALL * nx
            W = D_BALL * ny
            H = D_BALL * nz
            
            SA = 2 * (L*W + L*H + W*H)
            
            # Check if this configuration is an improvement
            if SA < SA0:
                # If we find a better box, record it
                if SA < min_found_sa:
                    min_found_sa = SA
                    best_config = (L, W, H, int(SA))

    # --- Output the result ---
    if best_config:
        L, W, H, SA = best_config
        # The problem asks to output each number in the final equation.
        # However, the format is a:b:c:d, not an equation.
        # We will print the formatted string as requested.
        print(f"{L}:{W}:{H}:{SA}")
    else:
        # If no better configuration is found, answer is 0.
        print(0)

solve()