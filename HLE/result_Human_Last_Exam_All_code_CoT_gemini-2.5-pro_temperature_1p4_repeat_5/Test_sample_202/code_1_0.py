import math

def find_shortest_distance(a, b, c):
    """
    Calculates the shortest distance between the maximum and minimum valued nodes
    in the described equilateral triangle setup.

    Args:
        a (float): The value at vertex A.
        b (float): The value at vertex B.
        c (float): The value at vertex C.
    """
    
    # The logic is based on the number of unique values among a, b, c.
    # We use a set to find the count of unique values.
    unique_values_count = len(set([a, b, c]))

    print(f"Given vertex values a={a}, b={b}, c={c}.")
    print(f"Number of unique vertex values is: {unique_values_count}")

    # Case 1: Three unique values. The max and min are at distinct vertices.
    # The distance is the side length of the triangle, which is 1.
    if unique_values_count == 3:
        distance = 1
        print("The maximum and minimum values are located at two different vertices.")
        print("The shortest distance is the side length of the triangle.")
        print(f"Final equation: distance = {distance}. The number in the equation is 1.")

    # Case 2: Two unique values. One extremum is a vertex, the other is the opposite side.
    # The distance is the altitude of the triangle, sqrt(3)/2.
    elif unique_values_count == 2:
        distance = math.sqrt(3) / 2
        print("One extreme value is at a vertex, and the other is along the opposite side.")
        print("The shortest distance is the altitude of the equilateral triangle.")
        # The instruction is to output each number in the final equation.
        # The equation is distance = sqrt(3)/2. The numbers are 3 and 2.
        print(f"Final equation: distance = sqrt(3)/2. The numbers in the equation are 3, 2.")
        print(f"The calculated distance is approximately {distance:.5f}")

    # Case 3: One unique value. All nodes have the same value.
    # The distance is 0.
    else: # unique_values_count == 1
        distance = 0
        print("All vertices have the same value, so all nodes have the same value.")
        print("The maximum and minimum values are equal.")
        print(f"Final equation: distance = {distance}. The number in the equation is 0.")

if __name__ == '__main__':
    # You can change these values to test different scenarios.
    # Case 1: a, b, c are all different.
    print("--- Case 1: Three unique values ---")
    find_shortest_distance(10, 20, 30)
    print("\n" + "="*40 + "\n")

    # Case 2: Two values are the same.
    print("--- Case 2: Two unique values ---")
    find_shortest_distance(15, 15, 5)
    print("\n" + "="*40 + "\n")

    # Case 3: All values are the same.
    print("--- Case 3: One unique value ---")
    find_shortest_distance(7, 7, 7)