import sys

def calculate_max_intersection(x, y):
    """
    Calculates the maximum intersection area for a given test square.

    The 4x4 square is decomposed into 16 1x1 squares (polygons).
    A test square S = [x, x+1] x [y, y+1] is checked against all 16 polygons.
    """
    if not (0 <= x <= 3 and 0 <= y <= 3):
        print("Error: x and y must be in the range [0, 3].")
        return

    max_area = 0.0
    
    # The 16 polygons are the unit squares [i, i+1] x [j, j+1]
    # for i, j in {0, 1, 2, 3}.
    for i in range(4):
        for j in range(4):
            # Polygon P_ij is the square [i, i+1] x [j, j+1]
            
            # Calculate the intersection of S and P_ij
            # S is [x, x+1] x [y, y+1]
            
            # Intersection in x-dimension
            x_overlap = max(0, min(x + 1, i + 1) - max(x, i))
            
            # Intersection in y-dimension
            y_overlap = max(0, min(y + 1, j + 1) - max(y, j))
            
            intersection_area = x_overlap * y_overlap
            
            if intersection_area > max_area:
                max_area = intersection_area
                
    print(f"For the test square at ({x}, {y}), the area is {max_area:.4f}")
    
# Let's test the worst-case scenario discussed above
x_test, y_test = 1.5, 1.5
print(f"Considering the simple grid decomposition into 16 unit squares.")
print(f"Let's test a square S centered at ({x_test+0.5}, {y_test+0.5}), with its bottom-left corner at ({x_test}, {y_test}).")
calculate_max_intersection(x_test, y_test)