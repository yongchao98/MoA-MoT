import math

def solve_intersection_volume():
    """
    Calculates the volume of the intersection of two shifted cubes.
    """
    # 1. Define the properties of the cubes and the shift.
    # Each cube has a side length of 1.
    side_length = 1.0
    
    # The shift is 1/2 along the main diagonal. We interpret this as
    # a shift vector of (1/2, 1/2, 1/2).
    shift_per_axis = 0.5
    
    print("This program calculates the volume of the intersection of two cubes.")
    print(f"Each cube has a side length of {side_length}.")
    print("-" * 30)

    # 2. Define the boundaries of the two cubes.
    # Let Cube 1 be at the origin, occupying the region [0,1] on each axis.
    c1_min = [0.0, 0.0, 0.0]
    c1_max = [side_length, side_length, side_length]
    print(f"Cube 1 is defined by the region where: 0 <= x,y,z <= {side_length}")

    # Cube 2 is shifted by `shift_per_axis` along each axis.
    c2_min = [c1_min[i] + shift_per_axis for i in range(3)]
    c2_max = [c1_max[i] + shift_per_axis for i in range(3)]
    print(f"Cube 2 is shifted by {shift_per_axis} along each axis to the region where: {c2_min[0]} <= x,y,z <= {c2_max[0]}")
    print("-" * 30)
    
    # 3. Calculate the boundaries of the intersection volume.
    # The intersection is also a cuboid. Its boundaries are determined by the
    # maximum of the minimum coordinates and the minimum of the maximum coordinates.
    intersect_min_x = max(c1_min[0], c2_min[0])
    intersect_min_y = max(c1_min[1], c2_min[1])
    intersect_min_z = max(c1_min[2], c2_min[2])

    intersect_max_x = min(c1_max[0], c2_max[0])
    intersect_max_y = min(c1_max[1], c2_max[1])
    intersect_max_z = min(c1_max[2], c2_max[2])
    
    print("The intersection forms a smaller cube with the following side lengths:")
    
    # 4. Calculate the side lengths of the intersection volume.
    side_x = intersect_max_x - intersect_min_x
    side_y = intersect_max_y - intersect_min_y
    side_z = intersect_max_z - intersect_min_z
    
    print(f"Side length (x-axis) = {intersect_max_x} - {intersect_min_x} = {side_x}")
    print(f"Side length (y-axis) = {intersect_max_y} - {intersect_min_y} = {side_y}")
    print(f"Side length (z-axis) = {intersect_max_z} - {intersect_min_z} = {side_z}")
    print("-" * 30)

    # 5. Calculate the final volume.
    volume = side_x * side_y * side_z

    print("The volume of this intersection is the product of its side lengths.")
    print(f"Final Equation: Volume = {side_x} * {side_y} * {side_z} = {volume}")

if __name__ == '__main__':
    solve_intersection_volume()