import math

def find_shortest_distance():
    """
    This script calculates the shortest distance between the node with the maximum
    number and the node with the minimum number on the described triangle grid.

    The solution is based on the insight that the value at any node is a linear
    interpolation of the values at the main vertices A, B, and C. This means the
    overall maximum and minimum values are always located at the boundaries of the
    large triangle (the vertices or edges).
    """
    try:
        # Prompt the user to enter the values for a, b, and c.
        input_str = input("Enter the three real numbers for vertices A, B, and C, separated by spaces: ")
        values = [float(x) for x in input_str.split()]
        if len(values) != 3:
            raise ValueError("Please provide exactly three numbers.")
        a, b, c = values[0], values[1], values[2]

    except (ValueError, IndexError) as e:
        print(f"Invalid input: {e}. Please run the script again and enter three numbers.")
        return

    # Count the number of unique values among a, b, and c.
    num_unique_values = len(set([a, b, c]))

    # Determine the shortest distance based on the number of unique values.
    if num_unique_values == 1:
        # Case 1: a = b = c. All nodes have the same value.
        print("\nAll three vertex values are the same (a = b = c).")
        print("Therefore, all nodes on the triangle have the same value.")
        print("The maximum and minimum values are equal, so the distance is 0.")
        final_distance = 0.0
        print(f"\nShortest Distance: {final_distance}")

    elif num_unique_values == 2:
        # Case 2: Two values are equal, one is different.
        # The maximum (or minimum) is along one edge, and the other extremum is at the opposite vertex.
        # The shortest distance is the altitude of the equilateral triangle.
        side = 1.0
        half_side = side / 2.0
        altitude_sq = side**2 - half_side**2
        altitude = math.sqrt(altitude_sq)
        
        print("\nExactly two of the vertex values are equal.")
        print("One extremum (max or min) lies on a vertex, and the other lies on the entire opposite edge.")
        print("The shortest distance is the altitude of the triangle with side length 1.")
        print("\nFinal Equation:")
        print(f"distance = sqrt(side^2 - (side/2)^2)")
        print(f"distance = sqrt({side}^2 - {half_side}^2) = sqrt({altitude_sq})")
        print(f"\nShortest Distance: {altitude}")

    else: # num_unique_values == 3
        # Case 3: All three values are different.
        # The maximum is at one vertex and the minimum is at another.
        # The distance is the side length of the triangle.
        side = 1.0
        print("\nAll three vertex values are distinct.")
        print("The maximum value is at one vertex and the minimum value is at another.")
        print("The shortest distance is the distance between these two vertices, which is the side length of the triangle.")
        final_distance = side
        print(f"\nShortest Distance: {final_distance}")

# Execute the main function
find_shortest_distance()