import math

def solve_triangle_grid_problem():
    """
    Calculates the largest number k of coordinate grid squares the triangle's perimeter can pass through.
    """
    # 1. Define the leg length of the right-angled isosceles triangle.
    L = 18

    # 2. Find the optimal orientation. This occurs when tan(theta) = 2.
    # We don't need the angle theta itself, but rather its sine and cosine.
    tan_theta = 2
    # sin(theta) = tan(theta) / sqrt(1 + tan(theta)^2)
    sin_theta = tan_theta / math.sqrt(1 + tan_theta**2)
    # cos(theta) = 1 / sqrt(1 + tan(theta)^2)
    cos_theta = 1 / math.sqrt(1 + tan_theta**2)

    # 3. Calculate the dimensions (Width and Height) of the triangle's bounding box
    # for this optimal orientation.
    # The vertices are at (0,0), (L*cos, L*sin), (-L*sin, L*cos).
    # The x-coordinates are 0, L*cos, -L*sin.
    # The y-coordinates are 0, L*sin, L*cos.
    # Width W = (L*cos) - (-L*sin) = L * (cos + sin)
    W = L * (cos_theta + sin_theta)
    # Height H = max(L*cos, L*sin) - 0. Since tan>1, sin>cos. So H = L*sin.
    H = L * sin_theta

    print(f"The leg length of the triangle is L = {L}.")
    print("The optimal orientation is found when tan(theta) = 2.")
    print(f"For this orientation, the bounding box width is W = {W:.4f}")
    print(f"And the bounding box height is H = {H:.4f}\n")

    # 4. The maximum integer spans for the coordinates are the floor of W and H.
    # By carefully placing the triangle (translation), we can achieve these spans.
    max_x_span = math.floor(W)
    max_y_span = math.floor(H)

    print(f"The maximum span of the integer parts of the x-coordinates is floor(W) = {max_x_span}.")
    print(f"The maximum span of the integer parts of the y-coordinates is floor(H) = {max_y_span}.\n")

    # 5. Calculate the maximum number of crossed squares, k.
    # The formula is k = 2 * (x_span) + 2 * (y_span).
    k = 2 * max_x_span + 2 * max_y_span

    print("The largest number of squares, k, is given by the formula:")
    print(f"k = 2 * {max_x_span} + 2 * {max_y_span}")
    print(f"k = {2 * max_x_span} + {2 * max_y_span}")
    print(f"k = {k}")


solve_triangle_grid_problem()
<<<80>>>