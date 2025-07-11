import math

def get_heesch_numbers():
    """
    Determines the Heesch number for the polygons in the image.

    Reasoning:
    1. All three polygons shown in the image are identical in shape.
    2. The polygon shape has point symmetry (it is unchanged after a 180-degree rotation about its center).
    3. The polygon can be classified as a centrally symmetric hexagon. Each of its six 'sides' is congruent to the opposite side after a 180-degree rotation.
       - The top-left jagged side is congruent to the bottom-right jagged side.
       - The top-right jagged side is congruent to the bottom-left jagged side.
       - The straight left side is congruent to the straight right side.
    4. A fundamental theorem in geometry states that any centrally symmetric hexagon tiles the plane.
    5. By definition, if a polygon can tile the plane, its Heesch number is considered to be infinity.
    6. Therefore, each of the three identical polygons has a Heesch number of infinity.
    """
    # Using the unicode character for infinity
    infinity_symbol = "\u221e"

    # The Heesch number is the same for all three identical polygons.
    heesch_number_1 = infinity_symbol
    heesch_number_2 = infinity_symbol
    heesch_number_3 = infinity_symbol

    print(f"{heesch_number_1}, {heesch_number_2}, {heesch_number_3}")

get_heesch_numbers()