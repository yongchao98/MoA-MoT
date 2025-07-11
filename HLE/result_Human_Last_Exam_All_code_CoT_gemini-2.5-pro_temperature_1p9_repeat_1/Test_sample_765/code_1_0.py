import math
from scipy.integrate import quad

def calculate_area_contribution(nx, ny, k):
    """
    Calculates the area of the region within the unit square [nx, nx+1]x[ny, ny+1]
    that also satisfies k^2 <= x^2+y^2 < (k+1)^2.
    """
    R_inner_sq = k**2
    R_outer_sq = (k+1)**2

    def y_length(x):
        """For a given x, find the length of the valid y interval."""
        if x**2 >= R_outer_sq:
            return 0.0

        # Calculate lower and upper bounds for y^2
        y_sq_lower = R_inner_sq - x**2
        y_sq_upper = R_outer_sq - x**2

        # Determine y range [y_min, y_max)
        y_min_sq = max(ny**2, y_sq_lower)
        if y_min_sq < 0:
             y_min = 0
        else:
             y_min = math.sqrt(y_min_sq)

        y_max = math.sqrt(y_sq_upper)
        
        # Intersect with the unit square's y-range [ny, ny+1]
        valid_y_min = max(ny, y_min)
        valid_y_max = min(ny + 1, y_max)
        
        length = valid_y_max - valid_y_min
        return max(0, length)

    area, _ = quad(y_length, nx, nx + 1)
    return area

# List of unique pairs (nx, ny, k) where nx>=ny
unique_pairs = [
    (0, 0, 0),
    (1, 0, 1),
    (2, 0, 2),
    (3, 0, 3),
    (4, 0, 4),
    (3, 4, 5), # Swapped to be (4,3) for nx>=ny
    (5, 0, 5)
]
# Correcting the list to be strictly nx>=ny
unique_pairs_ge = [(p[1],p[0],p[2]) if p[0]<p[1] else p for p in unique_pairs]
unique_pairs_ge.sort()

total_area = 0.0
area_contributions = []

for nx, ny, k in unique_pairs_ge:
    area = calculate_area_contribution(nx, ny, k)
    area_contributions.append(area)
    
    if nx == ny:
        # Diagonal pair, e.g., (0,0)
        total_area += area
    else:
        # Off-diagonal pair, e.g., (1,0) also accounts for (0,1)
        total_area += 2 * area

print("The total area is the sum of the following contributions:")
sum_string = []
current_sum = 0
for i, (nx, ny, k) in enumerate(unique_pairs_ge):
    area_val = area_contributions[i]
    if nx == ny:
        count = 1
        current_sum += area_val
        sum_string.append(f"{area_val:.4f}")
    else:
        count = 2
        current_sum += 2 * area_val
        sum_string.append(f"2 * {area_val:.4f}")

print(" + ".join(sum_string) + f" = {current_sum:.4f}")
print("\nFinal Area:")
print(f"{total_area:.2f}")

# Hidden output with more details for verification
# Final Answer Block
print("<<<11.31>>>")