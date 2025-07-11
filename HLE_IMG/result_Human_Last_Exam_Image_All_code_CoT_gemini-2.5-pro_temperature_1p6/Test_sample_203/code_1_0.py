import math

def solve_heesch_numbers():
    """
    This function analyzes the three polygons and states their Heesch numbers.
    """
    
    # Analysis of the three polygons
    p1_analysis = (
        "Polygon 1 (Top):\n"
        "This polygon is centrally symmetric. Each protrusion on one side has a corresponding, "
        "complementary indentation on the opposite side. Such shapes can tile the entire plane.\n"
        "Therefore, its Heesch number is infinity."
    )
    
    p2_analysis = (
        "\nPolygon 2 (Middle):\n"
        "This polygon has three 'protrusion' features and one complementary 'indentation' feature. "
        "This allows a full first layer (a corona) to be formed around a central tile, so H >= 1.\n"
        "However, the imbalance of 3 protrusions to 1 indentation means the boundary of this first "
        "corona is too 'spiky' to be surrounded by a second layer.\n"
        "Therefore, its Heesch number is 1."
    )
    
    p3_analysis = (
        "\nPolygon 3 (Bottom):\n"
        "This polygon has several protrusions but lacks the corresponding variety of indentations needed "
        "to fill them. For example, the 'W'-shaped protrusions on the top and bottom cannot be filled "
        "by any feature on another copy of the tile.\n"
        "Thus, the polygon cannot be completely surrounded even once.\n"
        "Therefore, its Heesch number is 0."
    )
    
    # The Heesch numbers for the three polygons in order (top, middle, bottom)
    h1 = float('inf')
    h2 = 1
    h3 = 0
    
    # Print the analysis
    print(p1_analysis)
    print(p2_analysis)
    print(p3_analysis)
    
    # Format the final answer string
    # Using 'infinity' for clear representation of the mathematical concept.
    answer_string = f"infinity, {h2}, {h3}"
    
    print("\n" + "="*40)
    print("The Heesch numbers for the three polygons in order are:")
    print(answer_string)
    print("="*40)

solve_heesch_numbers()
