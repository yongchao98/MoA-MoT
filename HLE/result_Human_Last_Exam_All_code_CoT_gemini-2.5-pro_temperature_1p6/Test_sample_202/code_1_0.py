import math

def solve_triangle_distance():
    """
    Calculates the shortest distance between the maximum and minimum numbered nodes
    in the described equilateral triangle.
    """
    # --- User-configurable values ---
    # Values at the vertices A, B, and C
    a = 10.0
    b = 50.0
    c = 10.0

    # Number of divisions of each side (this does not affect the final distance,
    # but is included for completeness based on the problem description).
    n = 10
    # --- End of user-configurable values ---

    # The vertices of the large equilateral triangle are A, B, C.
    # The side length of the triangle ABC is 1.
    side_length = 1

    print("Step 1: Understand the values at the nodes.")
    print("Based on the problem's rules, the value at any node is a linear interpolation of the vertex values a, b, and c.")
    print("The formula for a node P with integer coordinates (i, j, k) where i+j+k=n is:")
    print(f"V(P) = (a*i + b*j + c*k) / n")
    print(f"V(P) = ({a}*i + {b}*j + {c}*k) / {n}\n")


    print("Step 2: Find the maximum and minimum values.")
    max_val = max(a, b, c)
    min_val = min(a, b, c)
    print(f"The maximum possible value in the grid is max(a, b, c) = {max_val}.")
    print(f"The minimum possible value in the grid is min(a, b, c) = {min_val}.\n")

    print("Step 3: Determine the locations of the max/min values.")
    print("These extreme values occur at the main vertices A, B, or C of the triangle.")
    # Find which vertices have the max/min values
    vertices = {'A': a, 'B': b, 'C': c}
    max_locations = [v for v, val in vertices.items() if val == max_val]
    min_locations = [v for v, val in vertices.items() if val == min_val]
    print(f"The point(s) with the maximum value are at vertex/vertices: {', '.join(max_locations)}.")
    print(f"The point(s) with the minimum value are at vertex/vertices: {', '.join(min_locations)}.\n")


    print("Step 4: Calculate the shortest distance.")
    # If max and min values are the same, a, b, and c must be equal.
    # This means all nodes have the same value.
    if max_val == min_val:
        distance = 0
        print("Since the max and min values are equal, the points can be the same.")
        print("Shortest distance = 0")
    else:
        # If the values are different, the max and min points are at different vertices.
        # The distance between any two distinct vertices of the triangle is the side length.
        distance = side_length
        print("Since the max and min values are different, the points are at distinct vertices.")
        print(f"The distance between any two distinct vertices is the side length: {side_length}")

    print("\n---")
    print(f"Final Answer: The shortest distance is {distance}")
    # This is for the final answer tag format, assuming the code above ran.
    # We use this hidden way to prepare the final output string.
    globals()['final_answer_for_tag'] = distance

# Run the solver
solve_triangle_distance()
# The required final answer tag, populated by the function call.
# Using a specific example (a=10, b=50, c=10), max is 50, min is 10. They are different.
# So the distance is 1.
# final_answer_for_tag = 1
<<<1>>>