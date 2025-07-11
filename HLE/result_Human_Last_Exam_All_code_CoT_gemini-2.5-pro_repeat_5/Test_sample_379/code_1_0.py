import math

def solve_cube_locus_problem():
    """
    Calculates the length of a locus on a cube's surface and expresses it as a percentage.
    
    The problem asks to find the length of the locus of points C on the surface of a cube
    with side length r, where each point in C is at a surface distance r from a point P,
    which is the midpoint of an edge. The final answer should be the length of C divided
    by 2*pi*r, expressed as a whole number percentage.
    """
    
    # For calculation purposes, we can set r=1, as it will cancel out in the final ratio.
    r = 1.0
    pi = math.pi

    print("Step 1: Calculate the length of the locus on the Front face.")
    # The locus on the front face is an arc of a circle with radius r.
    # The angle subtended by this arc is pi/3 radians.
    len_front = (pi / 3) * r
    print(f"The arc on the Front face has length (pi/3) * r = {len_front:.4f}\n")

    print("Step 2: Calculate the length of the locus on the Right face.")
    # By symmetry, the arc on the right face has the same length.
    len_right = (pi / 3) * r
    print(f"The arc on the Right face has length (pi/3) * r = {len_right:.4f}\n")

    print("Step 3: Calculate the length of the locus on the Top face.")
    # The locus on the top face consists of two separate arcs. Each arc has length (pi/3)*r.
    len_top = 2 * (pi / 3) * r
    print(f"The locus on the Top face has total length (2*pi/3) * r = {len_top:.4f}\n")
    
    print("Step 4: Calculate the length of the locus on the Bottom face.")
    # By symmetry, the locus on the bottom face is the same as the top face.
    len_bottom = 2 * (pi / 3) * r
    print(f"The locus on the Bottom face has total length (2*pi/3) * r = {len_bottom:.4f}\n")

    print("Step 5: Calculate the total length of the locus C.")
    # The total length is the sum of the lengths on the four contributing faces.
    # Other faces (Left, Back) are too far for their main area to be part of the locus.
    total_length_C = len_front + len_right + len_top + len_bottom
    print(f"Total length C = {len_front:.4f} + {len_right:.4f} + {len_top:.4f} + {len_bottom:.4f} = {total_length_C:.4f}")
    print(f"In terms of pi and r, C = (pi/3 + pi/3 + 2*pi/3 + 2*pi/3)*r = (6*pi/3)*r = 2*pi*r\n")

    print("Step 6: Calculate the final value as a whole number percentage.")
    # The problem asks for (Length of C) / (2 * pi * r)
    denominator = 2 * pi * r
    ratio = total_length_C / denominator
    percentage = ratio * 100
    
    print("The final equation is: (Total Length C) / (2 * pi * r) * 100")
    print(f"Plugging in the numbers:")
    print(f"({total_length_C:.4f}) / ({denominator:.4f}) * 100 = {ratio:.2f} * 100 = {percentage:.0f}%")
    
    final_answer = int(percentage)
    print(f"\nThe final answer as a whole number percentage is {final_answer}.")
    
solve_cube_locus_problem()