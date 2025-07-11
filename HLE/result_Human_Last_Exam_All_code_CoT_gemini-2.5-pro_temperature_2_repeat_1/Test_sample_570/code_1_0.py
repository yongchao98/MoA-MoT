import math

def calculate_minimal_area():
    """
    This function calculates the minimal area of the convex domain K.
    The minimal domain is the rhombus (diamond shape) defined by |x| + |y| <= 1.
    """
    
    print("The convex domain with minimal area is the rhombus defined by |x| + |y| <= 1.")
    print("Its vertices are (1, 0), (0, 1), (-1, 0), and (0, -1).")
    
    # The diagonals of the rhombus connect opposite vertices.
    # The horizontal diagonal goes from (-1, 0) to (1, 0).
    # The vertical diagonal goes from (0, -1) to (0, 1).
    d1 = 1 - (-1)
    d2 = 1 - (-1)
    
    print(f"The length of the horizontal diagonal (d1) is {d1}.")
    print(f"The length of the vertical diagonal (d2) is {d2}.")
    
    # The area of a rhombus is calculated as (d1 * d2) / 2.
    area = (d1 * d2) / 2
    
    print("The area of a rhombus is given by the formula: (d1 * d2) / 2")
    print(f"So, the minimal area = ({d1} * {d2}) / 2 = {int(area)}")

calculate_minimal_area()
