import math

def get_heesch_numbers():
    """
    Determines the Heesch number for the polygons in the image.

    The polygon shown in the figure is a centrally symmetric hexagon.
    A key theorem in the study of tessellations states that any centrally
    symmetric hexagon can tile the plane (cover it without gaps or overlaps)
    using only translational movements.

    The Heesch number of a shape is the maximum number of layers of identical
    shapes that can be placed around it. If a shape can tile the entire plane,
    it can be surrounded by an infinite number of layers. By definition, its
    Heesch number is infinity.

    Since all three polygons in the image are identical, they all have the
    same Heesch number.
    """
    # The Heesch number for a shape that tiles the plane is infinity.
    heesch_number_1 = "infinity"
    heesch_number_2 = "infinity"
    heesch_number_3 = "infinity"

    # The problem asks for the answers separated by commas.
    print(f"{heesch_number_1}, {heesch_number_2}, {heesch_number_3}")

get_heesch_numbers()