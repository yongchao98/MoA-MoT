def solve_heesch_numbers():
    """
    This function determines and prints the Heesch numbers for the three polygons shown in the figure.

    The analysis is as follows:
    1.  The top polygon has complementary opposite sides, which allows it to tile the plane.
        Therefore, its Heesch number is infinite.
    2.  The middle polygon does not have fully complementary opposite sides. It is a known
        example of a shape with a Heesch number of 1.
    3.  The bottom polygon, similarly, does not tile the plane and is also a known
        example of a shape with a Heesch number of 1.
    """
    # Heesch numbers for the three polygons in order (top, middle, bottom)
    heesch_numbers = ["infinity", 1, 1]

    # The problem defines the Heesch number for a tiling polygon as infinity.
    polygon1_heesch = heesch_numbers[0]
    # The second and third polygons are known examples of shapes with Heesch number 1.
    polygon2_heesch = heesch_numbers[1]
    polygon3_heesch = heesch_numbers[2]
    
    # The question asks for the answers in order separated by commas.
    # The final equation/answer is the sequence of these numbers.
    print(f"{polygon1_heesch}, {polygon2_heesch}, {polygon3_heesch}")

solve_heesch_numbers()