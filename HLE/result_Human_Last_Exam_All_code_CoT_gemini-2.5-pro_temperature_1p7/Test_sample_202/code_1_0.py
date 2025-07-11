import math

def find_shortest_distance():
    # Define the given values for the vertices A, B, and C.
    # You can change these to test different scenarios.
    # Note: n, the number of divisions, does not affect the final distance.
    a = 10
    b = 10
    c = 5

    print(f"Analyzing the triangle with vertex values a={a}, b={b}, c={c}.")
    print("-" * 30)

    # --- Step 1 & 2: Determine node values and find overall max/min ---
    # The value at any node is a weighted average of a, b, and c.
    # This is derived from the rhombus rule, which implies a linear relationship.
    # Value(P) = u*a + v*b + w*c, where (u,v,w) are barycentric coordinates.
    # Therefore, the maximum value in the grid is max(a,b,c) and the minimum is min(a,b,c).
    
    values = sorted([a, b, c])
    min_val = values[0]
    mid_val = values[1]
    max_val = values[2]

    print(f"The minimum value across all nodes is: {min_val}")
    print(f"The maximum value across all nodes is: {max_val}")
    print("-" * 30)

    # --- Step 3 & 4: Locate extremal nodes and calculate distance ---
    # The shortest distance depends on how many of the values a,b,c are equal.

    # Case 1: All values are distinct (e.g., a > b > c).
    # The min is at one vertex, the max is at another vertex.
    if min_val < mid_val < max_val:
        print("All three vertex values are distinct.")
        print("The maximum value is at one vertex and the minimum is at another.")
        distance = 1.0
        print("The final equation for the distance is simply the side length.")
        print("The number in this equation is 1.")
        print(f"\nShortest distance = {distance}")

    # Case 2: All values are equal (a = b = c).
    # The max and min values are the same and occur everywhere.
    elif min_val == max_val:
        print("All three vertex values are equal.")
        print("The maximum and minimum values are the same.")
        distance = 0.0
        print("The final equation for the distance is 0.")
        print("The number in this equation is 0.")
        print(f"\nShortest distance = {distance}")
        
    # Case 3: Exactly two values are equal (e.g., a = b > c or a > b = c).
    # One extremum (min or max) is at a single vertex.
    # The other extremum is located along the entire opposite edge.
    # The shortest distance is the altitude of the equilateral triangle.
    else:
        print("Exactly two of the vertex values are equal.")
        print("One extremum (max or min) is at a vertex, the other is on the opposite edge.")
        # Altitude of equilateral triangle with side 1 is sqrt(1^2 - (1/2)^2) = sqrt(3)/2.
        distance = math.sqrt(3) / 2
        print("The final equation for the distance is sqrt(3) / 2.")
        print("The numbers in this equation are 3 and 2.")
        print(f"\nShortest distance = {distance}")

    return distance

# Execute the function to find and print the result.
final_answer = find_shortest_distance()
print(f"<<<{final_answer}>>>")