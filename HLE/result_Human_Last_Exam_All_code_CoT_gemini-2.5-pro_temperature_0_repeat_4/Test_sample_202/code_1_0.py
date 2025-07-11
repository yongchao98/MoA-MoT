import math

def solve_triangle_distance():
    """
    Calculates the shortest distance between the max and min value nodes
    in a described equilateral triangle grid.
    """
    # In an equilateral triangle ABC with side length 1, the vertices A, B, C
    # have values a, b, c.
    # Please provide the values for a, b, and c.
    # Example values are provided below.
    a = 5
    b = 10
    c = 5

    print("Step 1: Define the values at the vertices A, B, and C.")
    print(f"Value at A (a): {a}")
    print(f"Value at B (b): {b}")
    print(f"Value at C (c): {c}\n")

    print("Step 2: Analyze the location of the maximum and minimum values.")
    print("The value at any node is a linear interpolation of a, b, and c.")
    print("Therefore, the absolute maximum and minimum values must be located at the vertices A, B, or C, or on the edges connecting them if vertex values are equal.\n")

    values = [a, b, c]
    max_val = max(values)
    min_val = min(values)

    print(f"Step 3: Determine the shortest distance based on the values.")
    # Case 1: All vertex values are the same.
    if max_val == min_val:
        distance = 0.0
        print("Analysis: All vertex values are equal (a=b=c).")
        print("Every node in the triangle has the same value.")
        print("The shortest distance between a point with the max value and a point with the min value is 0.")
        print("\nFinal Equation: distance = 0")

    else:
        count_max = values.count(max_val)
        count_min = values.count(min_val)

        # Case 2: Two vertices share the max value, or two share the min value.
        if count_max >= 2 or count_min >= 2:
            side = 1.0
            distance = math.sqrt(side**2 - (side/2)**2)
            print("Analysis: Two vertices share an extreme value (e.g., a=b>c or a>b=c).")
            print("One extreme value is found along an entire edge of the triangle.")
            print("The other extreme value is at the opposite vertex.")
            print("The shortest distance is the altitude of the triangle.")
            print(f"\nFinal Equation: distance = sqrt(side^2 - (side/2)^2)")
            print(f"distance = sqrt({side}^2 - ({side/2})^2) = sqrt({side**2 - (side/2)**2}) = {distance}")

        # Case 3: All vertex values are distinct.
        else:
            distance = 1.0
            print("Analysis: The vertex values a, b, and c are all distinct.")
            print("The maximum value is at one vertex, and the minimum is at another.")
            print("The shortest distance is the side length of the triangle.")
            print("\nFinal Equation: distance = 1")

    print("\n" + "="*40)
    print(f"The shortest distance is: {distance}")
    print("="*40)

# Execute the function
solve_triangle_distance()