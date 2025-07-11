import sys

def solve_temporal_fixed_points():
    """
    Finds and sums all positive integer temporal fixed points less than 100
    based on the provided rules.
    """
    # The temporal horizon constant.
    h = 48
    # The upper limit for the fixed points.
    limit = 100

    # A list to store the found fixed points.
    fixed_points = []

    # A value 'k' is a temporal fixed point if |k+k| = |k-k|.
    # With the temporal horizon h(T) = 48, this is interpreted as:
    # |2k| mod 48 = |0| mod 48
    # Since we are looking for positive integers k, this simplifies to (2 * k) % 48 == 0.
    for k in range(1, limit):
        if (2 * k) % h == 0:
            fixed_points.append(k)

    # Calculate the sum of the found fixed points.
    total_sum = sum(fixed_points)

    # Build the equation string to display the final calculation.
    # The list of numbers is converted to strings, joined by ' + ',
    # and the sum is appended.
    # For example: "24 + 48 + 72 + 96 = 240"
    sum_equation = " + ".join(map(str, fixed_points)) + " = " + str(total_sum)

    # Print the final result in the required format.
    print("The temporal fixed points are the positive integers k < 100 where (2 * k) is a multiple of 48.")
    print(f"Found fixed points: {fixed_points}")
    print("The sum of these fixed points is calculated as follows:")
    print(sum_equation)

solve_temporal_fixed_points()
<<<240>>>