import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-sided regular polygon
    to the 2n-sided regular polygon it was constructed from.

    The problem states: "start with a 2n sided regular hexagon and extend
    alternate edges until they intersect to form a n sided regular polygon".
    For a hexagon, 2n=6, so n=3. The resulting polygon is a triangle.
    """

    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    # The general formula for the ratio is (1 + sec(pi/n)) / 2
    # Let's calculate this step-by-step for the given n.

    print(f"Calculating the area ratio for n = {n}:")
    print(f"The starting polygon has 2n = {2*n} sides.")
    print(f"The constructed polygon has n = {n} sides.")
    print("-" * 30)

    # Step 1: Calculate theta = pi / n in radians
    theta = math.pi / n
    print(f"1. The angle (pi/n) is: pi / {n} = {theta:.4f} radians.")

    # Step 2: Calculate cos(theta)
    cos_theta = math.cos(theta)
    print(f"2. The cosine of this angle is: cos({theta:.4f}) = {cos_theta:.4f}.")

    # Step 3: Calculate sec(theta)
    sec_theta = 1 / cos_theta
    print(f"3. The secant is the reciprocal of the cosine: sec({theta:.4f}) = {sec_theta:.4f}.")

    # Step 4: Apply the formula (1 + sec(theta)) / 2
    # We output each number in the final equation as requested.
    numerator_val = 1 + sec_theta
    denominator_val = 2
    ratio = numerator_val / denominator_val

    print(f"4. The final ratio is (1 + sec(pi/n)) / 2.")
    print(f"   Using the values: ({1} + {sec_theta:.4f}) / {denominator_val} = {ratio:.4f}")
    print("-" * 30)
    print(f"Result: The area is {ratio} times larger.")


# We use n=3 as in the problem's example (hexagon -> triangle)
calculate_area_ratio(3)