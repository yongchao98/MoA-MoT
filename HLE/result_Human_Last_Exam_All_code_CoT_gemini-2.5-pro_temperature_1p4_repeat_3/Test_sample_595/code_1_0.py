import math

def calculate_squares(p1, p2):
    """Calculates the number of squares a line segment passes through."""
    # Using the formula for non-axis-aligned segments that don't cross grid vertices
    # Number of squares = N_v + N_h + 1
    # N_v = |floor(x2) - floor(x1)|
    # N_h = |floor(y2) - floor(y1)|
    
    x1, y1 = p1
    x2, y2 = p2
    
    # Handle axis-aligned cases separately
    if abs(x1 - x2) < 1e-9: # Vertical
        return abs(math.floor(y2) - math.floor(y1)) + 1
    if abs(y1 - y2) < 1e-9: # Horizontal
        return abs(math.floor(x2) - math.floor(x1)) + 1
        
    nv = abs(math.floor(x2) - math.floor(x1))
    nh = abs(math.floor(y2) - math.floor(y1))
    
    return nv + nh + 1

# --- Triangle Properties ---
leg = 18
hyp = 18 * math.sqrt(2)

# Optimal angle found to be arctan(1/2)
angle = math.atan(1/2)
cos_a = math.cos(angle)
sin_a = math.sin(angle)

# --- Optimal Placement ---
# Choose a small epsilon to offset from integer coordinates
epsilon = 0.01

# Vertices
B = (epsilon, epsilon)

# Vector BA
dx_ba = leg * cos_a
dy_ba = leg * sin_a
A = (B[0] + dx_ba, B[1] + dy_ba)

# Vector BC
dx_bc = -leg * sin_a
dy_bc = leg * cos_a
C = (B[0] + dx_bc, B[1] + dy_bc)

# Calculate squares for each side
s_ba = calculate_squares(B, A)
s_bc = calculate_squares(B, C)
s_ac = calculate_squares(A, C)

# Total unique squares
# Overlaps occur at the 3 vertices, assuming they are in distinct squares
# B is in (0,0), A is in (16,8), C is in (-9,16) -> distinct
k = s_ba + s_bc + s_ac - 3

print(f"Optimal Orientation Angle: arctan(1/2) ~ {math.degrees(angle):.2f} degrees")
print(f"Vertex B: ({B[0]:.4f}, {B[1]:.4f})")
print(f"Vertex A: ({A[0]:.4f}, {A[1]:.4f})")
print(f"Vertex C: ({C[0]:.4f}, {C[1]:.4f})")
print("-" * 30)
print(f"Squares crossed by side BA: {s_ba}")
print(f"Squares crossed by side BC: {s_bc}")
print(f"Squares crossed by side AC: {s_ac}")
print(f"Total squares (sum): {s_ba + s_bc + s_ac}")
print("Subtracting 3 for the squares containing the vertices (overlaps)...")
print(f"k = {s_ba} + {s_bc} + {s_ac} - 3 = {k}")
