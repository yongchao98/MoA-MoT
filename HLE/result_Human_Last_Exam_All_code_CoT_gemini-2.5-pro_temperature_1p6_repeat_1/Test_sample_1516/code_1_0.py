import math

def solve_parliament_design():
    """
    Calculates the maximum integer value for the parameter K in the parliament design.

    The key constraint is that the head of any seated observer must be at a height
    greater than or equal to the head of the standing speaker. This ensures the line
    of sight does not slope downwards.

    H_seated(j) >= H_standing(0)
    (r_j^2 / K) + 1 >= (r_0^2 / K) + 1.5

    This must hold for all rows j=1 to 12. The most restrictive case is for the
    first row of observers (j=1), which gives the lowest upper bound for K.

    K <= 2 * (r_1^2 - r_0^2)
    """

    # Initial parameters
    r_0 = 3.0  # Radial distance of the first row (m)
    row_depth = 1.5  # Depth of each row (m)
    
    # Calculate the radial distance of the second row (first observer row, j=1)
    r_1 = r_0 + 1 * row_depth

    # Calculate the maximum allowed value for K
    k_max_float = 2 * (r_1**2 - r_0**2)
    
    # K must be an integer
    k_max_int = math.floor(k_max_float)

    # Print the explanation and the equation with its values
    print("The problem is to find the maximum integer value for K.")
    print("The governing constraint is that every seated member's head height must be >= the standing speaker's head height.")
    print("This leads to the inequality: K <= 2 * (r_1^2 - r_0^2)")
    print("\nCalculating the values for the equation:")
    print(f"The multiplier is: 2")
    print(f"The value for r_1 (radial distance of the first observer row) is: {r_1}")
    print(f"The value for r_0 (radial distance of the speaker's row) is: {r_0}")
    
    print("\nSubstituting the values into the equation:")
    print(f"K <= 2 * ({r_1}^2 - {r_0}^2)")
    print(f"K <= 2 * ({r_1**2} - {r_0**2})")
    print(f"K <= 2 * ({r_1**2 - r_0**2})")
    print(f"K <= {k_max_float}")
    
    print(f"\nSince K must be an integer, the maximum value is the floor of {k_max_float}.")
    print(f"The maximum integer value K can take is: {k_max_int}")

solve_parliament_design()
<<<22>>>