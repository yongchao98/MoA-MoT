import math

# --- Define Constants ---

# Dimensions of the material
MATERIAL_W = 140
MATERIAL_H = 110

# Dimensions and character counts for the artifacts
SQUARE_SIDE = 10
SQUARE_CHARS = 4

CIRCLE_RADIUS = 20
CIRCLE_DIAMETER = CIRCLE_RADIUS * 2
CIRCLE_CHARS = 999

def calculate_optimal_layout(material_w, material_h):
    """
    Calculates the optimal number of circles and squares to maximize characters.
    This function implements a greedy approach: it prioritizes circles due to
    their higher character-per-area value, then fills the rest with squares.
    """
    # --- Step 1: Maximize the number of Circles (M) ---
    # We fit as many 40x40 bounding boxes for circles as possible.
    cols_circles = math.floor(material_w / CIRCLE_DIAMETER)
    rows_circles = math.floor(material_h / CIRCLE_DIAMETER)
    num_circles = cols_circles * rows_circles

    # --- Step 2: Calculate remaining area for Squares (N) ---
    # The circles occupy a rectangular block at one corner of the material.
    used_width = cols_circles * CIRCLE_DIAMETER
    used_height = rows_circles * CIRCLE_DIAMETER

    # The remaining area is composed of two rectangular strips.
    # Strip 1:
    rem_strip1_w = material_w - used_width
    rem_strip1_h = material_h
    squares_in_strip1 = math.floor(rem_strip1_w / SQUARE_SIDE) * math.floor(rem_strip1_h / SQUARE_SIDE)

    # Strip 2:
    rem_strip2_w = used_width
    rem_strip2_h = material_h - used_height
    squares_in_strip2 = math.floor(rem_strip2_w / SQUARE_SIDE) * math.floor(rem_strip2_h / SQUARE_SIDE)

    num_squares = squares_in_strip1 + squares_in_strip2
    
    return num_squares, num_circles

# --- Main Calculation ---

# We check both orientations of the material to ensure we find the absolute maximum.
# Orientation 1: 140x110
n1, m1 = calculate_optimal_layout(MATERIAL_W, MATERIAL_H)
k1 = (SQUARE_CHARS * n1) + (CIRCLE_CHARS * m1)

# Orientation 2: 110x140 (rotated)
n2, m2 = calculate_optimal_layout(MATERIAL_H, MATERIAL_W)
k2 = (SQUARE_CHARS * n2) + (CIRCLE_CHARS * m2)

# Select the best result (in this case, they are identical)
if k1 > k2:
    N = n1
    M = m1
    K = k1
else:
    N = n2
    M = m2
    K = k2

# --- Output the result ---
print("To maximize the number of engraved characters, the workers should produce:")
print(f"- {N} squares for the name.")
print(f"- {M} circles for the bio.")
print("\nThis yields a maximum total of K characters, calculated as follows:")
print(f"K = ({SQUARE_CHARS} * {N}) + ({CIRCLE_CHARS} * {M})")
print(f"K = {SQUARE_CHARS * N} + {CIRCLE_CHARS * M}")
print(f"K = {K}")
print(f"\nThe final answer in the format N:M:K is: {N}:{M}:{K}")

# The final answer for the system
final_answer_string = f"<<<{N}:{M}:{K}>>>"