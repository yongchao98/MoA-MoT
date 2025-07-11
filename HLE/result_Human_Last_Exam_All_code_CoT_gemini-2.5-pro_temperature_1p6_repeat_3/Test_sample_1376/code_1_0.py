def solve_tomb_artifacts():
    """
    This script solves the problem of maximizing the number of engraved Chinese characters
    by finding the optimal number of squares (N) and circles (M) that can be
    cut from a given rectangular piece of meteorite material.
    """
    # Step 1: Define problem parameters
    plate_w = 140  # cm
    plate_h = 110  # cm
    square_dim = 10  # cm
    circle_dim = 20 * 2  # cm, diameter of the circle
    chars_per_square = 4

    # Step 2: Determine characters per circle based on encoding requirements
    # A circle holds 999 symbols. The bio has 7*7*9 = 441 unique characters.
    # To encode 441 states with 8 unique Bagua symbols, we need 3 symbols
    # per character (since 8^2=64 < 441 and 8^3=512 > 441).
    chars_per_circle = 999 // 3

    # Step 3 & 4: Execute the packing strategy: maximize circles first, then fill with squares.
    # Calculate the maximum number of circles (M) that can fit on the plate.
    m_count_w = plate_w // circle_dim  # floor(140 / 40) = 3
    m_count_h = plate_h // circle_dim  # floor(110 / 40) = 2
    M = m_count_w * m_count_h

    # Calculate the area occupied by the circles
    circles_area_w = m_count_w * circle_dim # 3 * 40 = 120 cm
    circles_area_h = m_count_h * circle_dim # 2 * 40 = 80 cm

    # Calculate how many squares (N) fit in the remaining L-shaped area.
    # The L-shape is split into two rectangular sections.
    # Section 1:
    rem_section1_w = plate_w - circles_area_w # 140 - 120 = 20
    rem_section1_h = circles_area_h          # 80
    squares_in_section1 = (rem_section1_w // square_dim) * (rem_section1_h // square_dim)

    # Section 2:
    rem_section2_w = plate_w                  # 140
    rem_section2_h = plate_h - circles_area_h # 110 - 80 = 30
    squares_in_section2 = (rem_section2_w // square_dim) * (rem_section2_h // square_dim)
    
    N = squares_in_section1 + squares_in_section2

    # Step 5: Calculate the maximum total characters (K)
    K = (chars_per_square * N) + (chars_per_circle * M)
    
    # Print the final results as requested
    print("To maximize the number of engraved characters, the optimal combination is:")
    print(f"N (number of squares) = {N}")
    print(f"M (number of circles) = {M}")
    print("")
    print("The final result in the format N:M:K is:")
    print(f"{N}:{M}:{K}")
    print("")
    print("The final calculation for the total number of characters (K) is:")
    print(f"{chars_per_square} * {N} + {chars_per_circle} * {M} = {K}")

solve_tomb_artifacts()