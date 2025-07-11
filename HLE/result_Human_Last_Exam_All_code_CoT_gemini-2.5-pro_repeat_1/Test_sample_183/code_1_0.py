import math

def solve():
    """
    This function calculates the number of integers n in the range [10^18, 10^18 + 10000]
    that can be expressed in the form n = x^3 + 2y^3 + 4z^3 - 6xyz for some integers x, y, z.

    The logic, as outlined in the plan, focuses on the case x = 10^6 and y=0 or z=0.
    """

    # Case 1: y = 0. We need 0 <= 4z^3 <= 10000, so 0 <= z^3 <= 2500.
    # We count the number of positive integers z satisfying this.
    z_max_cubed = 2500
    # Calculate the maximum integer z such that z^3 <= 2500.
    z_max = int(math.pow(z_max_cubed, 1/3))
    count_z_positive = z_max

    # Case 2: z = 0. We need 0 <= 2y^3 <= 10000, so 0 <= y^3 <= 5000.
    # We count the number of positive integers y satisfying this.
    y_max_cubed = 5000
    # Calculate the maximum integer y such that y^3 <= 5000.
    y_max = int(math.pow(y_max_cubed, 1/3))
    count_y_positive = y_max

    # Case 3: y = 0 and z = 0. This gives n = 10^18.
    count_base = 1

    # The total number of unique integers is the sum of counts from the three disjoint cases:
    # 1. y > 0, z = 0
    # 2. z > 0, y = 0
    # 3. y = 0, z = 0
    total_count = count_y_positive + count_z_positive + count_base

    # Print the result as an equation, as requested.
    print(f"{count_z_positive} + {count_y_positive} + {count_base} = {total_count}")

solve()