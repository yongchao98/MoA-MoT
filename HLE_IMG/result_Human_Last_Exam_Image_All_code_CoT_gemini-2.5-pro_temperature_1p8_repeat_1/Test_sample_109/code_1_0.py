import math

def solve_target_problem():
    """
    Calculates the expected number of hits on yellow circles in 10000 shots.
    """
    # Number of shots
    num_shots = 10000

    # Properties of the yellow circles
    num_yellow_circles = 6
    r_yellow = 1  # in cm

    # Properties of the white circles, deduced from the problem statement
    # R^2 + (R-1)^2 = (R+1)^2 => R^2 - 4R = 0 => R = 4
    R_white = 4 # in cm

    # To satisfy the "all coordinates are integers" constraint, the layout must be a grid.
    # The image's staggered layout is geometrically inconsistent with this constraint.
    # For a grid layout, the horizontal shift dx=0 and the vertical shift dy=2*R=8.
    
    # Calculate the total target dimensions H and W
    # The vertical distance between centers dy = H/2 - R => 8 = H/2 - 4 => H = 24
    H = 24 # in cm
    # The width is deduced to be the space of 4 large circles side-by-side.
    W = 4 * (2 * R_white) # in cm

    # Calculate total area of the target
    total_area = W * H

    # Calculate total area of the yellow circles
    total_yellow_area = num_yellow_circles * math.pi * r_yellow**2

    # Calculate the probability of hitting a yellow circle
    probability = total_yellow_area / total_area

    # Calculate the expected number of hits
    expected_hits = num_shots * probability

    # Print the equation step by step
    print("The calculation for the expected number of hits is as follows:")
    print(f"Expected Hits = Total Shots * (Total Yellow Area / Total Target Area)")
    print(f"Equation: Expected Hits = {num_shots} * ({num_yellow_circles} * pi * {r_yellow}^2) / ({W} * {H})")
    print(f"               = {num_shots} * ({num_yellow_circles} * pi * {r_yellow**2}) / ({total_area})")
    print(f"               = {num_shots * num_yellow_circles} * pi / {total_area}")
    
    # Simplify the fraction for the final expression
    numerator = num_shots * num_yellow_circles
    denominator = total_area
    common_divisor = math.gcd(numerator, denominator)
    
    print(f"               = ({numerator // common_divisor} * pi) / {denominator // common_divisor}")
    print("\nFinal Answer:")
    print(f"The expected number of hits is {expected_hits:.1f}")

solve_target_problem()
<<<245.4>>>