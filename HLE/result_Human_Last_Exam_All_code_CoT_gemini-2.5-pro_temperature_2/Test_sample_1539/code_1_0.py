import math

def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters.
    """

    # --- Step 1: Define Properties ---
    material_w = 140  # cm
    material_h = 110  # cm

    # Properties of the "Book of Heaven" circular plate
    circle_radius = 20
    circle_diameter = circle_radius * 2
    circle_symbols = 9999

    # Properties of the "Book of Earth" squared plate
    square_side = 10
    square_symbols = 360

    # --- Step 2: Packing Strategy ---
    # The value density of circles is much higher (9999 symbols / (40*40cm area) = 6.25 symbols/cm^2)
    # compared to squares (360 symbols / (10*10cm area) = 3.6 symbols/cm^2).
    # Therefore, the strategy is to maximize the number of circles first.

    # --- Step 3: Calculate Maximum Circular Plates (N) ---
    # We treat circles as 40x40 squares for packing purposes (guillotine cuts).
    # Consider both orientations of the material.
    
    # Orientation 1: 140x110
    circles_fit_w1 = math.floor(material_w / circle_diameter) # floor(140/40) = 3
    circles_fit_h1 = math.floor(material_h / circle_diameter) # floor(110/40) = 2
    N1 = circles_fit_w1 * circles_fit_h1 # 3 * 2 = 6

    # Orientation 2: 110x140
    circles_fit_w2 = math.floor(material_h / circle_diameter) # floor(110/40) = 2
    circles_fit_h2 = math.floor(material_w / circle_diameter) # floor(140/40) = 3
    N2 = circles_fit_w2 * circles_fit_h2 # 2 * 3 = 6
    
    # The maximum number of circles is 6.
    N = N1
    
    # --- Step 4: Calculate Maximum Square Plates (M) in Remaining Area ---
    # Let's use the first orientation. The 6 circles occupy a block of:
    circles_block_w = circles_fit_w1 * circle_diameter  # 3 * 40 = 120 cm
    circles_block_h = circles_fit_h1 * circle_diameter  # 2 * 40 = 80 cm

    # The remaining area consists of two rectangular strips.
    # Strip 1:
    rem_strip1_w = material_w - circles_block_w  # 140 - 120 = 20 cm
    rem_strip1_h = material_h                    # 110 cm
    squares_in_strip1 = math.floor(rem_strip1_w / square_side) * math.floor(rem_strip1_h / square_side)
    # (20/10) * (110/10) = 2 * 11 = 22

    # Strip 2:
    rem_strip2_w = circles_block_w              # 120 cm
    rem_strip2_h = material_h - circles_block_h  # 110 - 80 = 30 cm
    squares_in_strip2 = math.floor(rem_strip2_w / square_side) * math.floor(rem_strip2_h / square_side)
    # (120/10) * (30/10) = 12 * 3 = 36

    M = squares_in_strip1 + squares_in_strip2
    
    # --- Step 5: Calculate Total Symbols (K) ---
    K_from_circles = N * circle_symbols
    K_from_squares = M * square_symbols
    K = K_from_circles + K_from_squares

    # --- Step 6: Format the Output ---
    print("This solution maximizes the number of engraved symbols by prioritizing the more valuable circular plates.")
    print(f"\nNumber of circular plates (N): {N}")
    print(f"Number of squared plates (M): {M}")
    print("\nCalculating the maximal number of engraved characters (K):")
    print(f"K = ({N} circular plates * {circle_symbols} symbols) + ({M} squared plates * {square_symbols} symbols)")
    print(f"K = {K_from_circles} + {K_from_squares}")
    print(f"K = {K}")
    print(f"\nThe final answer is in the format N:M:K")
    print(f"<<<{N}:{M}:{K}>>>")

solve_laozi_engraving()