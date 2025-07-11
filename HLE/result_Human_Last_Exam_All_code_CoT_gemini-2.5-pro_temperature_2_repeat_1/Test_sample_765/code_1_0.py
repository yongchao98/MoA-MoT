import math
from scipy.integrate import quad

def get_intersection_area(R_squared, x_min, y_min, x_max, y_max):
    """
    Calculates the area of intersection of a circle x^2+y^2 < R_squared 
    and a rectangle [x_min, x_max] x [y_min, y_max], assuming the rectangle 
    is in the first quadrant.
    """
    if R_squared <= 0:
        return 0
    R = math.sqrt(R_squared)

    def integrand(x):
        # We are only interested in x within the circle's radius
        if x >= R:
            return 0
        
        y_circle = math.sqrt(R_squared - x*x)
        
        # Calculate the height of the intersecting slice at position x
        height = max(0, min(y_circle, y_max) - y_min)
        return height

    # The integration needs to be over the rectangle's x-range
    # capped by the circle's radius R.
    integration_domain_max = min(x_max, R)
    if x_min >= integration_domain_max:
        return 0

    area, _ = quad(integrand, x_min, integration_domain_max)
    return area

# Main calculation logic
total_area = 0
terms = []
term_labels = []

# Find all integer pairs (n,m) in 0<=n,m<6 where sqrt(n^2+m^2) is an integer
valid_pairs = []
for n in range(6):
    for m in range(6):
        k_val = math.sqrt(n**2 + m**2)
        if k_val == int(k_val):
            valid_pairs.append({'n': n, 'm': m, 'k': int(k_val)})

# Sum up the areas from all valid squares
for pair in valid_pairs:
    n, m, k = pair['n'], pair['m'], pair['k']

    # The region is defined by k^2 <= a^2+b^2 < (k+1)^2.
    # The condition a^2+b^2 >= k^2 is always satisfied for a point in [n, n+1)x[m,m+1).
    # So we only need to calculate the area for a^2+b^2 < (k+1)^2.
    
    R_squared = (k + 1)**2
    
    # Define unit square boundaries
    x_min, x_max = n, n + 1
    y_min, y_max = m, m + 1
    
    area_contribution = get_intersection_area(R_squared, x_min, y_min, x_max, y_max)

    total_area += area_contribution
    terms.append(area_contribution)
    term_labels.append(f"A({n},{m})")

# Print the final equation
equation_parts = [f"{area:.4f}" for area in terms]
print("The total area is the sum of the areas from 13 unit squares:")
print(" + ".join(equation_parts) + f" = {total_area:.4f}")
print("\nFinal Area:")
# Express the answer to two decimals
print(f"{total_area:.2f}")
