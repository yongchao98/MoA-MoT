import math

def solve_symmetry():
    """
    This function determines the rotational symmetry of the provided tiling.
    """
    
    # Step 1: Analyze the tiling to find a center of rotation.
    # By observing the tiling, we can choose the center of a dark blue hexagon
    # as a center of rotation.
    
    # Step 2: Count the number of times the pattern repeats around this center.
    # Around the center of a blue hexagon, we can see the pattern is built from
    # 6 identical surrounding motifs. Each motif consists of a yellow square
    # and an orange rhombus, along with the adjacent blue hexagons.
    # Rotating the entire image by a certain angle maps each motif onto the next one.
    # The number of motifs gives the order of rotational symmetry.
    order_n = 6
    
    # Step 3: Calculate the smallest angle of rotation that leaves the pattern unchanged.
    # For an n-fold rotational symmetry, the angle is 360/n degrees.
    full_circle = 360
    angle_of_rotation = full_circle / order_n
    
    # Step 4: Print the reasoning and the final result.
    print(f"The tiling exhibits rotational symmetry.")
    print(f"By choosing the center of a blue hexagon as the rotation axis, we observe {order_n} identical repeating units around it.")
    print(f"This indicates a {order_n}-fold rotational symmetry.")
    print("The angle for this symmetry operation is calculated as:")
    print(f"{full_circle} degrees / {order_n} = {int(angle_of_rotation)} degrees.")
    print(f"\nThe rotational symmetry of the tiling is {order_n}.")

solve_symmetry()
