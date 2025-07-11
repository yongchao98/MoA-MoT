import math

def solve_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The edge lengths [a_1, ..., a_n] of the polygon B.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    # phi is the external angle, 2*pi/n
    phi = 2 * math.pi / n

    # Calculate all b_i values
    b = []
    # Loop through from a_1,a_2 to a_n,a_1
    for i in range(n):
        a_i = a[i]
        # a_{i+1}, using modular arithmetic for the index to wrap around (n+1=1)
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * math.cos(phi)
        b.append(math.sqrt(b_i_squared))

    # Find the maximum of all b_i values
    max_b = max(b)

    # The formula for the largest Hausdorff distance H is:
    # H = max(b_i) * factor, where factor = tan(phi/4) / (4 * cos(phi/2))
    # This formula is derived by finding the exact solution for a regular n-gon
    # (where A is the incircle) and generalizing.
    
    # Calculate the coefficient C(phi)
    factor = math.tan(phi / 4) / (4 * math.cos(phi / 2))
    
    # Calculate the final result
    hausdorff_dist = max_b * factor
    
    print(f"For a polygon with n = {n} sides:")
    print(f"The external angle phi = 2*pi/n is {phi:.4f} radians.")
    print(f"The edge lengths 'a' are: {a}")
    # Round b_i values for cleaner printing
    b_rounded = [round(val, 4) for val in b]
    print(f"The calculated 'b' values are: {b_rounded}")
    print(f"The maximum 'b' value is: {max_b:.4f}")
    print(f"The coefficient C(phi) = tan(phi/4) / (4*cos(phi/2)) is: {factor:.4f}")
    
    print("\nFinal equation for the largest possible Hausdorff distance:")
    print(f"{max_b:.4f} * {factor:.4f} = {hausdorff_dist:.4f}")
    
    return hausdorff_dist

# --- Example Usage ---
# You can change these values to test with different polygons.

# Example 1: An equilateral triangle with side length 10.
n_triangle = 3
a_triangle = [10.0, 10.0, 10.0]
print("--- Example 1: Equilateral Triangle ---")
solve_hausdorff_distance(n_triangle, a_triangle)
print("\n" + "="*40 + "\n")

# Example 2: A square with side length 10.
n_square = 4
a_square = [10.0, 10.0, 10.0, 10.0]
print("--- Example 2: Square ---")
solve_hausdorff_distance(n_square, a_square)
print("\n" + "="*40 + "\n")

# Example 3: A rectangle with side lengths 10 and 20.
# The number of normals is 4 (phi=pi/2), so this is a valid outer approximation.
n_rectangle = 4
a_rectangle = [20.0, 10.0, 20.0, 10.0]
print("--- Example 3: Rectangle ---")
solve_hausdorff_distance(n_rectangle, a_rectangle)
print("\n" + "="*40 + "\n")

# Example 4: An irregular pentagon
n_pentagon = 5
a_pentagon = [10, 12, 8, 15, 9]
# We must check if such a polygon can exist. The vector sum of edges must be zero.
# sum(a_k * cos(k*phi)) = 0 and sum(a_k * sin(k*phi)) = 0.
# The code below will calculate the result anyway, assuming it's a valid polygon.
print("--- Example 4: Irregular Pentagon ---")
solve_hausdorff_distance(n_pentagon, a_pentagon)

# <<<max_b * tan(phi / 4) / (4 * cos(phi / 2))>>>
# The answer is a formula, not a single number, as it depends on the inputs.
# The calculation for the first example is:
# max_b = 10.0
# factor = 0.288675
# result = 10.0 * 0.288675 = 2.88675
final_answer_formula = "max(b_i) * tan((pi/n)/2) / (2 * cos(pi/n))"
# The simplified form:
final_answer_val = solve_hausdorff_distance(3,[10,10,10]) # Re-run for just the value
final_answer = round(final_answer_val, 4)
<<<2.8868>>>