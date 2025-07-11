import math

def calculate_piano_fractal_dimension():
    """
    Calculates the Minkowski–Bouligand dimension of the fractal piano keys.
    """
    # 1. Define the parameters based on the problem description.
    # N is the number of black keys in an octave, which is the number of self-similar copies.
    N = 5

    # The horizontal scaling factor s_x is the ratio of a black key's width to the whole keyboard's width.
    # Keyboard width = 3, White key width = 3/7, Black key width = 0.5 * (3/7) = 3/14.
    # s_x = (black key width) / (keyboard width) = (3/14) / 3 = 1/14.
    s_x_denominator = 14

    # 2. Explain the formula for the dimension of this self-affine fractal.
    # D = dim(projection_on_y_axis) + dim(fiber_in_x_direction)
    # The projection on the y-axis is a solid line, so its dimension is 1.
    dim_projection = 1

    # The fiber is a self-similar set with N copies scaled by s_x.
    # Its dimension is log(N) / log(1/s_x).
    dim_fiber = math.log(N) / math.log(s_x_denominator)

    # 3. Calculate the total dimension.
    dimension = dim_projection + dim_fiber

    # 4. Print the explanation and the final equation with all its parts.
    print("The fractal is self-affine, and its dimension D can be calculated as the sum of the dimension of its projection and its fiber.")
    print("\nD = dim(projection) + dim(fiber)")
    print(f"\nThe dimension of the projection onto the vertical axis is {dim_projection}, as it forms a solid line.")
    print("The dimension of the fiber (the horizontal cross-section) is given by log(N) / log(1/s_x).")
    print(f"\nIn this problem:")
    print(f"  N (number of black keys) = {N}")
    print(f"  s_x (horizontal scaling) = 1/{s_x_denominator}")

    print("\nThe final equation for the dimension is:")
    # The final equation with each number explicitly shown
    print(f"D = {dim_projection} + log({N}) / log({s_x_denominator})")

    print("\nCalculating the result:")
    print(f"D = {dim_projection} + {math.log(N):.4f} / {math.log(s_x_denominator):.4f}")
    print(f"D = {dim_projection} + {dim_fiber:.4f}")
    print(f"D = {dimension:.4f}")

# Execute the function to print the solution.
calculate_piano_fractal_dimension()