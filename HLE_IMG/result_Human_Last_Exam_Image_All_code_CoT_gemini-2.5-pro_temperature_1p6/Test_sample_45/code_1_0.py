import math

def solve_symmetry():
    """
    This function determines and explains the rotational symmetry of the tiling.
    """
    # The order of rotational symmetry is determined by counting the number of times
    # a motif repeats around a central point.
    # By observing the star-like patterns, we can see a central heptagon (7-sided polygon)
    # with 7 arms radiating from it.
    order_of_symmetry = 7
    
    # The angle of rotation is 360 degrees divided by the order of symmetry.
    angle = 360 / order_of_symmetry
    
    print(f"The tiling has a rotational symmetry of order {order_of_symmetry}.")
    print("This can be seen by focusing on one of the central star-shaped motifs.")
    print("The motif repeats itself 7 times around a central point.")
    print("Therefore, rotating the tiling by the following angle leaves it unchanged:")
    print(f"360 / {order_of_symmetry} = {angle:.2f} degrees")

solve_symmetry()