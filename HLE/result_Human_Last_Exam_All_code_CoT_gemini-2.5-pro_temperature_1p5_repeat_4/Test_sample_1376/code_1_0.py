import math

def solve_emperor_tomb_problem():
    """
    Calculates the optimal number of squares and circles to maximize the
    total number of engraved Chinese characters from a single meteorite sheet.
    """
    # --- Step 1: Define constants and problem parameters ---

    # A set of 7x7x9 unique Chinese characters for the bio
    unique_chars_for_bio = 7 * 7 * 9  # 441

    # Bagua has 8 symbols. Find how many symbols are needed per character.
    # We need to solve for s in 8^s >= 441.
    # 8^2 = 64 (not enough)
    # 8^3 = 512 (sufficient)
    symbols_per_char = 3

    # Define the number of characters per artifact
    chars_per_square = 4  # The 4 characters of "Qin Shi Huang Di"
    symbols_per_circle = 999
    # Calculate how many full characters fit on a circle
    chars_per_circle = math.floor(symbols_per_circle / symbols_per_char)  # 333

    # Define artifact dimensions for cutting (using a square bounding box for the circle)
    square_dim = 10  # 10x10 cm
    circle_dim = 40  # 20cm radius -> 40cm diameter

    # Define the two possible orientations of the meteorite sheet
    sheet_orientations = [(140, 110), (110, 140)]

    # --- Step 2: Systematically search for the optimal layout ---

    best_N = 0
    best_M = 0
    max_K = 0

    # Function to find the best layout for a given sheet orientation
    def find_best_layout(sheet_W, sheet_H):
        local_best_N = 0
        local_best_M = 0
        local_max_K = 0

        # Maximum number of circles that can fit along the width and height
        max_circles_W = math.floor(sheet_W / circle_dim)
        max_circles_H = math.floor(sheet_H / circle_dim)

        # Iterate through all possible grid arrangements of circles
        for m_w in range(max_circles_W + 1):
            for m_h in range(max_circles_H + 1):
                M = m_w * m_h

                # Calculate the area used by the grid of circles
                width_used_by_circles = m_w * circle_dim
                height_used_by_circles = m_h * circle_dim

                # The remaining area is L-shaped. We can split it into two rectangles to
                # calculate the number of squares that can be cut.
                
                # Area 1: The vertical part of the 'L'
                rem_width1 = sheet_W - width_used_by_circles
                rem_height1 = sheet_H
                n1 = math.floor(rem_width1 / square_dim) * math.floor(rem_height1 / square_dim)

                # Area 2: The horizontal part of the 'L' (below the circles)
                rem_width2 = width_used_by_circles
                rem_height2 = sheet_H - height_used_by_circles
                n2 = math.floor(rem_width2 / square_dim) * math.floor(rem_height2 / square_dim)

                N = n1 + n2
                K = (chars_per_square * N) + (chars_per_circle * M)

                if K > local_max_K:
                    local_max_K = K
                    local_best_N = N
                    local_best_M = M
        
        return local_best_N, local_best_M, local_max_K

    # Check both sheet orientations and find the overall best result
    for W, H in sheet_orientations:
        N, M, K = find_best_layout(W, H)
        if K > max_K:
            max_K = K
            best_N = N
            best_M = M
            
    # --- Step 3: Print the final results and the calculation ---
    print(f"To maximize the number of characters, the optimal combination is N = {best_N} squares and M = {best_M} circles.")
    print("\nThe maximal number of Chinese characters (K) is calculated as follows:")
    print(f"K = (Characters per Square * N) + (Characters per Circle * M)")
    print(f"K = ({chars_per_square} * {best_N}) + ({chars_per_circle} * {best_M})")
    print(f"K = {chars_per_square * best_N} + {chars_per_circle * best_M}")
    print(f"K = {max_K}")
    print(f"\nThe final answer in the format N:M:K is: {best_N}:{best_M}:{max_K}")

# Execute the function to solve the problem
solve_emperor_tomb_problem()