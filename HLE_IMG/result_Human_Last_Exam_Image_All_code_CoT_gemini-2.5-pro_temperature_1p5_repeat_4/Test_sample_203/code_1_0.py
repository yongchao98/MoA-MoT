import math

def solve_heesch_number():
    """
    This function determines the Heesch number for the polygons in the image.
    
    1.  **Analyze the Shape:** The polygon shown is an irregular hexagon built on a triangular grid.
        - The top side has an inward V-shaped 'notch'.
        - The bottom side has an outward V-shaped 'tab'.
        - The top-left side has a 'notch'.
        - The bottom-right side has a 'tab'.
        - The top-right and bottom-left sides are straight.
    
    2.  **Check for Tiling:** The 'tab' on each side is shaped to perfectly fit into the 'notch' on the opposite side of an adjacent, identical polygon. This is a key feature of many shapes that tile the plane (tessellate). A hexagon with opposite sides modified in this complementary way will tile the plane.
    
    3.  **Define Heesch Number:** The Heesch number H(P) of a polygon P is the maximum number of layers of non-overlapping copies of P that can surround the original polygon. If a polygon can tile the plane, it can be surrounded by an infinite number of layers. By definition, its Heesch number is infinity.
    
    4.  **Conclusion:** Since the polygon tiles the plane, its Heesch number is infinity. The image shows three identical polygons, so they all share this property.
    """
    
    # Unicode for the infinity symbol
    infinity_symbol = "\u221e"
    
    # Heesch number for each of the three identical polygons
    heesch_number_1 = infinity_symbol
    heesch_number_2 = infinity_symbol
    heesch_number_3 = infinity_symbol
    
    # The final answer requires the three numbers separated by commas.
    print(f"{heesch_number_1}, {heesch_number_2}, {heesch_number_3}")

solve_heesch_number()