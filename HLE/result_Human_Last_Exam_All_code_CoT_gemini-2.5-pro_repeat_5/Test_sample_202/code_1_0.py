def find_shortest_distance():
    """
    Calculates the shortest distance between the node with the maximum value
    and the node with the minimum value in the described triangle.
    """
    # We will use example values for a, b, c, and n to demonstrate the solution.
    a = 10.0
    b = 20.0
    c = 30.0
    n = 100

    print("--- Problem Setup ---")
    print(f"Values at vertices: a={a} (at A), b={b} (at B), c={c} (at C)")
    print(f"Side divisions n = {n}")
    print(f"Triangle side length = 1")
    print("-" * 23)

    # Step 1: Define the equation for the value at any node (i, j).
    # As derived in the plan, the value is a linear function: f(i, j) = K1*i + K2*j + K0.
    # We calculate the coefficients K0, K1, K2.
    K0 = a
    K1 = (b - a) / n
    K2 = (c - a) / n

    # Step 2: Output the numbers in the final equation, as requested.
    print("\nThe value at a node (i, j) is given by the equation: f(i, j) = K1*i + K2*j + K0")
    print("The numbers in the final equation are:")
    print(f"K0 (value at A) = {K0}")
    print(f"K1 (slope along AB) = ({b} - {a}) / {n} = {K1}")
    print(f"K2 (slope along AC) = ({c} - {a}) / {n} = {K2}")
    
    # Step 3: Find the maximum and minimum values on the grid.
    # These must occur at the vertices A, B, or C.
    max_val = max(a, b, c)
    min_val = min(a, b, c)

    print(f"\nThe maximum value on any node is max(a,b,c) = {max_val}")
    print(f"The minimum value on any node is min(a,b,c) = {min_val}")

    # Step 4: Determine the shortest distance based on the max and min values.
    if max_val == min_val:
        # If a=b=c, all nodes have the same value. The distance is 0.
        shortest_distance = 0.0
        print("\nSince the maximum and minimum values are equal, the distance is 0.")
    else:
        # If the values are different, the max and min points are distinct
        # vertices of the main triangle. The distance between them is 1.
        shortest_distance = 1.0
        print("\nSince the max and min values are different, the points are distinct vertices.")
        print("The distance between any two distinct vertices of a triangle with side length 1 is 1.")

    print("\n--- Final Answer ---")
    print(f"The shortest distance is: {shortest_distance}")

# Execute the function
find_shortest_distance()