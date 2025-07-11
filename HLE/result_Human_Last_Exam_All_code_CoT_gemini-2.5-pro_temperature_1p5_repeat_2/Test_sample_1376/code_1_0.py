import math

# --- Define Constants ---

# Dimensions of the meteorite material rectangle
RECT_W = 140
RECT_H = 110

# Dimensions of the artifacts
SQUARE_SIDE = 10
# A circle with 20cm radius requires a 40x40cm square to be cut from
CIRCLE_BOX_SIDE = 40 

# Number of characters per artifact
CHARS_PER_SQUARE = 4
CHARS_PER_CIRCLE = 999

# --- Step 1: Maximize the number of circles (M) ---

# Calculate how many 40x40 squares fit into the 140x110 rectangle
# We consider both orientations, although in this case they yield the same result.
m_option1 = math.floor(RECT_W / CIRCLE_BOX_SIDE) * math.floor(RECT_H / CIRCLE_BOX_SIDE)
m_option2 = math.floor(RECT_H / CIRCLE_BOX_SIDE) * math.floor(RECT_W / CIRCLE_BOX_SIDE)
M = max(m_option1, m_option2)

# --- Step 2: Calculate the number of squares (N) in the remaining area ---

# We arrange the M circles in a dense block to maximize the remaining contiguous area.
# The block of circles will occupy an area of:
# (floor(RECT_W / CIRCLE_BOX_SIDE) * CIRCLE_BOX_SIDE) x (floor(RECT_H / CIRCLE_BOX_SIDE) * CIRCLE_BOX_SIDE)
w_used_by_circles = math.floor(RECT_W / CIRCLE_BOX_SIDE) * CIRCLE_BOX_SIDE
h_used_by_circles = math.floor(RECT_H / CIRCLE_BOX_SIDE) * CIRCLE_BOX_SIDE

# This leaves two rectangular areas.
# Area 1: The strip along the full width
rem_area1_w = RECT_W
rem_area1_h = RECT_H - h_used_by_circles

# Area 2: The strip along the side of the used height
rem_area2_w = RECT_W - w_used_by_circles
rem_area2_h = h_used_by_circles

# Calculate how many 10x10 squares fit in each remaining rectangular area
n1 = math.floor(rem_area1_w / SQUARE_SIDE) * math.floor(rem_area1_h / SQUARE_SIDE)
n2 = math.floor(rem_area2_w / SQUARE_SIDE) * math.floor(rem_area2_h / SQUARE_SIDE)

# The total number of small squares is the sum from both areas
N = n1 + n2

# --- Step 3: Calculate the total number of characters (K) ---

K = (N * CHARS_PER_SQUARE) + (M * CHARS_PER_CIRCLE)

# --- Print the results ---
print(f"To maximize the number of characters, the workers should produce:")
print(f"Number of 10x10cm squares (N): {N}")
print(f"Number of 20cm radius circles (M): {M}")
print("\nThe maximal number of Chinese characters (K) is calculated as:")
print(f"K = {N} * {CHARS_PER_SQUARE} + {M} * {CHARS_PER_CIRCLE} = {K}")
print(f"\nThe final answer in the format N:M:K is: {N}:{M}:{K}")

# Final answer in the required format
final_answer = f"<<<{N}:{M}:{K}>>>"
print(final_answer)
