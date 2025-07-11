import math

def find_shortest_distance():
    """
    Calculates the shortest distance between the maximum and minimum valued nodes
    on the described triangular grid.
    """
    try:
        # Prompt the user for the values at vertices A, B, and C.
        a = float(input("Enter the number on vertex A (a): "))
        b = float(input("Enter the number on vertex B (b): "))
        c = float(input("Enter the number on vertex C (c): "))
    except ValueError:
        print("Invalid input. Please enter valid numbers.")
        return

    # The maximum and minimum values on the grid will be the max and min of a, b, c.
    # We need to find the number of unique values among a, b, and c to determine
    # the location of the max/min points.
    unique_values = set([a, b, c])
    num_unique = len(unique_values)

    print("-" * 30)

    if num_unique == 1:
        # Case 1: All values are the same (a=b=c).
        # The max and min values are the same, located everywhere.
        print(f"The values a, b, c are all equal to {a}.")
        print("The shortest distance between the point with the maximum number and the point with the minimum number is 0.")
        print("The final equation for the distance is d = 0.")
        print("The number in this equation is 0.")

    elif num_unique == 2:
        # Case 2: Two values are the same, one is different (e.g., a=b != c).
        # One extremum is at a vertex, the other is along the opposite edge.
        # The distance is the altitude of the equilateral triangle.
        distance = math.sqrt(3) / 2
        print(f"The values a, b, c ({a}, {b}, {c}) have two unique values.")
        print(f"The shortest distance is the altitude of the triangle, which is sqrt(3)/2 ≈ {distance:.6f}.")
        print("The final equation for the distance is d = sqrt(3)/2.")
        print("The numbers in this equation are 3 and 2.")

    elif num_unique == 3:
        # Case 3: All values are different.
        # The max is at one vertex and the min is at another.
        # The distance is the side length of the triangle.
        distance = 1.0
        print(f"The values a, b, c ({a}, {b}, {c}) are all distinct.")
        print(f"The shortest distance is the side length of the triangle, which is {distance}.")
        print("The final equation for the distance is d = 1.")
        print("The number in this equation is 1.")

if __name__ == "__main__":
    find_shortest_distance()
