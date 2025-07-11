import math

def find_shortest_distance():
    """
    Calculates the shortest distance between the nodes with the maximum and
    minimum values on the described triangle grid.

    The solution depends on the relationship between the values a, b, and c
    at the main vertices A, B, and C.
    """

    # --- INPUTS ---
    # You can change these values to test different scenarios.
    # Case 1 example: a=10, b=5, c=1 (3 distinct values)
    # Case 2 example: a=10, b=10, c=1 (2 distinct values)
    # Case 3 example: a=5, b=5, c=5 (1 distinct value)
    a = 10.0
    b = 10.0
    c = 1.0
    
    print(f"Analyzing for vertex values: a={a}, b={b}, c={c}")

    # Determine the number of unique values among a, b, c
    unique_values = len(set([a, b, c]))

    # --- CALCULATION ---
    if unique_values == 1:
        # Case 3: a = b = c
        # The value is the same everywhere. Max and Min points are the same.
        distance = 0.0
        print("All vertex values are equal.")
        print("The set of max value points and min value points are identical.")
        print(f"Final shortest distance = {distance}")

    elif unique_values == 2:
        # Case 2: Two values are equal (e.g., a = b != c)
        # One extremum is on an edge, the other is at the opposite vertex.
        # The shortest distance is the altitude of the triangle.
        side = 1.0
        half_side = 0.5
        distance = math.sqrt(side**2 - half_side**2)
        print("Two vertex values are equal.")
        print("The shortest distance is the altitude of the triangle.")
        print(f"Equation: distance = sqrt(side^2 - (side/2)^2)")
        print(f"Numbers in equation: side={side}, half_side={half_side}")
        print(f"Final shortest distance = {distance}")

    else: # unique_values == 3
        # Case 1: a, b, c are all distinct
        # The max is at one vertex, the min is at another.
        # The shortest distance is the side length between two vertices.
        side = 1.0
        distance = side
        print("All vertex values are distinct.")
        print("The shortest distance is between two vertices of the triangle.")
        print(f"Equation: distance = side_length")
        print(f"Numbers in equation: side_length={side}")
        print(f"Final shortest distance = {distance}")

# Execute the function to see the result for the provided a,b,c values.
find_shortest_distance()