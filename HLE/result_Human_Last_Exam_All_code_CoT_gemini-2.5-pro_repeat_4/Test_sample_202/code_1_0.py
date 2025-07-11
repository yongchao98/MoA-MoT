import math

def solve_triangle_distance(a, b, c, n):
    """
    Solves for the shortest distance between the node with the maximum number
    and the node with the minimum number in the described triangle grid.

    The function explains the reasoning step-by-step.
    """

    print("--- Problem Analysis ---")
    print(f"We are given an equilateral triangle ABC with side length 1.")
    print(f"The side is divided into n={n} parts.")
    print(f"The values at the vertices are: f(A) = {a}, f(B) = {b}, f(C) = {c}.")
    print("\nThe crucial property is that for any small rhombus, the sum of values at opposite vertices is equal.")
    print("This implies that the value at any node is a linear interpolation of the vertex values a, b, and c.")
    print("Specifically, for a point P with barycentric coordinates (u, v, w) with respect to (A, B, C), its value is f(P) = u*a + v*b + w*c.")
    print("\n--- Finding Extrema ---")
    print("Because the value function is a linear combination (a weighted average) of a, b, and c, its maximum and minimum values over the entire triangle must be max(a, b, c) and min(a, b, c) respectively.")
    max_val = max(a, b, c)
    min_val = min(a, b, c)
    print(f"The maximum value on the grid is: {max_val}")
    print(f"The minimum value on the grid is: {min_val}")

    print("\n--- Locating Extrema and Calculating Distance ---")
    
    # Use a set to count the number of unique values among a, b, c
    unique_values = set([a, b, c])
    num_unique = len(unique_values)

    if num_unique == 1:
        print("Case 1: The values at all three vertices are equal (a=b=c).")
        print("This means every node in the triangle has the same value.")
        print("The points with maximum and minimum value are the same.")
        distance_str = "0"
        final_distance = 0
    elif num_unique == 3:
        print("Case 2: The values at the vertices a, b, and c are all distinct.")
        print("The maximum value occurs at one vertex, and the minimum at another.")
        print("The distance between any two distinct vertices of a unit equilateral triangle is 1.")
        distance_str = "1"
        final_distance = 1
    else: # num_unique == 2
        print("Case 3: Two vertex values are identical, and the third is different.")
        print("One extremum (max or min) occurs at a single vertex.")
        print("The other extremum occurs along the entire edge connecting the other two vertices.")
        print("The shortest distance is therefore the altitude from the single vertex to the opposite edge.")
        print("For a unit equilateral triangle, the altitude h is calculated as:")
        # Outputting the equation as requested
        print("h = sqrt(side^2 - (side/2)^2) = sqrt(1^2 - (1/2)^2)")
        final_distance = math.sqrt(1**2 - (0.5)**2)
        distance_str = f"sqrt(3)/2 ≈ {final_distance:.4f}"

    print("\n--- Final Answer ---")
    print(f"The shortest distance between the point with the maximum number and the point with the minimum number is: {distance_str}")
    return final_distance

# --- Execution with an example ---
# Let's choose values for a, b, c to demonstrate a case.
# n is part of the problem statement but not needed for the final calculation.
example_a = 7
example_b = 2
example_c = 7
example_n = 10

solve_triangle_distance(example_a, example_b, example_c, example_n)
