import math

def calculate_shortest_distance():
    """
    This function prompts the user for the values at the vertices of the triangle,
    calculates the shortest distance between the node with the maximum value and 
    the node with the minimum value, and prints the result.
    """
    # --- Get User Input ---
    while True:
        try:
            a_str = input("Enter the number on point A: ")
            a = float(a_str)
            b_str = input("Enter the number on point B: ")
            b = float(b_str)
            c_str = input("Enter the number on point C: ")
            c = float(c_str)
            break
        except ValueError:
            print("\nInvalid input. Please enter valid real numbers.\n")

    # --- Case Analysis ---
    values = [a, b, c]
    
    # Case 1: All values are equal.
    # The set of max points and min points are the same (all nodes).
    # The shortest distance is 0.
    if a == b and b == c:
        distance = 0.0
    else:
        max_val = max(values)
        min_val = min(values)
        
        count_max = values.count(max_val)
        count_min = values.count(min_val)
        
        # Case 2: Two vertices share the max value OR two share the min value.
        # This means the set of max points is an edge and the min point is a vertex
        # (or vice-versa). The shortest distance is the triangle's altitude.
        # For a side length of 1, altitude h = sqrt(1^2 - (1/2)^2) = sqrt(3)/2.
        if count_max == 2 or count_min == 2:
            distance = math.sqrt(3) / 2
        else:
            # Case 3: The values a, b, c are all distinct.
            # The max point is one vertex and the min point is another.
            # The shortest distance is the side length of the triangle.
            distance = 1.0
            
    # --- Print the Result ---
    print("\n-------------------------------------------------")
    print(f"For vertex values a={a}, b={b}, c={c}:")
    print(f"The shortest distance between the maximum and minimum value nodes is: {distance}")
    print("-------------------------------------------------")

# Execute the function
calculate_shortest_distance()