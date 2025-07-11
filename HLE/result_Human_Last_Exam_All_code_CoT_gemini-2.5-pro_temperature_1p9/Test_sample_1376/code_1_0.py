import math

def solve_maximization_problem():
    """
    Solves the artifact maximization problem for Emperor Qin.
    This function calculates the optimal number of squares (N) and circles (M)
    that can be produced from a rectangular meteorite material to maximize the
    total number of engraved characters (K).
    """

    # --- Problem Parameters ---

    # Material dimensions (cm)
    bin_width = 140
    bin_height = 110

    # Square artifact properties
    square_dim = 10
    square_value = 4  # characters

    # Circle artifact properties
    circle_radius = 20
    circle_bounding_box_dim = circle_radius * 2
    circle_value = 999  # characters (symbols)

    best_config = {
        "N": 0,
        "M": 0,
        "K": -1,
        "orientation": ""
    }

    # Consider both orientations of the material
    orientations = [
        (bin_width, bin_height),
        (bin_height, bin_width)
    ]

    for w, h in orientations:
        # --- Step 1: Maximize the number of high-value circles (M) ---
        # We pack them in a simple grid, which is optimal for identical items.
        
        # Number of circles that can fit along width and height
        m_fit_w = w // circle_bounding_box_dim
        m_fit_h = h // circle_bounding_box_dim
        M = m_fit_w * m_fit_h

        # Area occupied by the circles' bounding boxes
        occupied_w = m_fit_w * circle_bounding_box_dim
        occupied_h = m_fit_h * circle_bounding_box_dim

        # --- Step 2: Use leftover material for squares (N) ---
        # The leftover space is an L-shape. We can divide it into two rectangles
        # to calculate the number of squares that can be cut.
        
        # Rectangle 1: The vertical strip of leftover material
        rect1_w = w - occupied_w
        rect1_h = h
        n1 = (rect1_w // square_dim) * (rect1_h // square_dim)

        # Rectangle 2: The horizontal strip of leftover material (below the circles' area)
        rect2_w = occupied_w
        rect2_h = h - occupied_h
        n2 = (rect2_w // square_dim) * (rect2_h // square_dim)

        N = n1 + n2
        
        # --- Step 3: Calculate the total value (K) ---
        K = (N * square_value) + (M * circle_value)

        # Check if this configuration is better than the best one found so far
        if K > best_config["K"]:
            best_config["N"] = N
            best_config["M"] = M
            best_config["K"] = K
            best_config["orientation"] = f"{w}x{h}"

    # --- Step 4: Output the final result ---
    N = best_config['N']
    M = best_config['M']
    K = best_config['K']
    
    print(f"To maximize the number of engraved characters, the workers should produce:")
    print(f"- {N} squares for the official name (N)")
    print(f"- {M} circles for the bio (M)")
    print("\nThis yields a maximal total character count (K).")
    
    # Calculate parts of the equation for printing
    square_chars = N * square_value
    circle_chars = M * circle_value
    
    print("\nThe final calculation is:")
    print(f"K = (Characters from squares) + (Characters from circles)")
    print(f"K = ({square_value} * N) + ({circle_value} * M)")
    print(f"K = ({square_value} * {N}) + ({circle_value} * {M})")
    print(f"K = {square_chars} + {circle_chars}")
    print(f"K = {K}")

    print("\nThe answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")


solve_maximization_problem()
<<<58:6:6226>>>