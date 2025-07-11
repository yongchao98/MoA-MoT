import math

def find_shortest_distance():
    """
    Calculates the shortest distance between the node with the maximum number
    and the node with the minimum number in the described triangle grid.
    
    The user should modify the values of a, b, and c below.
    """
    # The numbers placed on points A, B, C are a, b, c respectively.
    # Please modify these values as needed.
    a = 10
    b = 20
    c = 10
    
    print(f"Given values: a = {a} (at vertex A), b = {b} (at vertex B), c = {c} (at vertex C)\n")

    # Case 1: All three values are equal.
    if a == b and b == c:
        print("Analysis: The values at all three vertices are equal.")
        print("This means all nodes in the triangle have the same value.")
        print("The set of maximum points and minimum points are the same.")
        print("\nFinal equation: distance = 0")
        final_distance = 0.0
        
    # Case 2: Exactly two of the three values are equal.
    elif a == b or b == c or a == c:
        print("Analysis: Exactly two vertices have the same value.")
        print("One extremum (max or min) is located at a single vertex.")
        print("The other extremum is located along the entire opposite side.")
        print("The shortest distance is the altitude of the equilateral triangle.")
        
        print("\nFinal equation: distance = sqrt(3) / 2")
        print("Numbers in the equation: 3, 2")
        final_distance = math.sqrt(3) / 2
        
    # Case 3: All three values are distinct.
    else:
        print("Analysis: The values at the three vertices are all distinct.")
        print("The maximum value is located at one vertex.")
        print("The minimum value is located at another vertex.")
        print("The shortest distance is the length of the side connecting these two vertices.")
        
        print("\nFinal equation: distance = 1")
        print("Numbers in the equation: 1")
        final_distance = 1.0

    print(f"\nResulting shortest distance: {final_distance}")

if __name__ == '__main__':
    find_shortest_distance()