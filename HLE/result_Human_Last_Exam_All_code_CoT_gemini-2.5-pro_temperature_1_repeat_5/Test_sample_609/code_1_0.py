import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-sided polygon to a 2n-sided polygon
    constructed by extending alternate edges of the 2n-gon.

    Args:
        n (int): The number of sides of the outer polygon.
    """
    if n <= 2:
        print("n must be greater than 2.")
        return

    # The general formula for the ratio of the areas is: (1 + sec(pi/n)) / 2
    # Let's calculate the components of this formula.

    # 1. Calculate theta = pi / n
    theta = math.pi / n

    # 2. Calculate sec(theta). math.sec doesn't exist, so we use 1/math.cos().
    cos_theta = math.cos(theta)
    if abs(cos_theta) < 1e-9:
        print(f"The formula is undefined for n={n} as cos(pi/n) is zero.")
        return
    sec_theta = 1 / cos_theta

    # 3. Calculate the final ratio using the formula.
    ratio = (1 + sec_theta) / 2

    # 4. Print the results in the required format.
    print(f"For a starting {2*n}-sided polygon, we construct a {n}-sided polygon.")
    print(f"The area of the new polygon is {ratio:.4f} times larger than the original polygon.")
    print("\nThe final equation with the calculated numbers is:")
    
    # We use formatting to make the output clear, matching the example's result of 1.5.
    print(f"(1 + {sec_theta:.1f}) / 2 = {ratio:.1f}")

# Run the calculation for the case given in the problem (a hexagon, 2n=6, so n=3).
calculate_area_ratio(3)
