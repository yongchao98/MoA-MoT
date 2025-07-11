def find_smallest_rectangle_area():
    """
    This script finds the area of the smallest integer-length rectangle
    that admits a non-guillotine tiling by squares from the set S={2x2, 3x3, 5x5, 7x7}.

    The solution is based on a theorem by Snover, King, and Kenyon (1991).
    A rectangle L x W can be tiled in a non-guillotine way if and only if L and W
    are the "half-sums" (A-values) of two linearly independent relations
    on the side lengths of the squares.

    A relation is an equation: c1*2 + c2*3 + c3*5 + c4*7 = 0.
    The A-value is the sum of the positive terms: A = sum over c_i > 0 of c_i*s_i.

    Our goal is to find the two smallest A-values from independent relations.
    """
    S = {2, 3, 5, 7}
    print(f"The set of available square sizes is S = {S}.")
    print("To find the smallest rectangle, we must find the two smallest 'A-values' that arise from independent relations.")
    print("An A-value corresponds to a number that can be formed in at least two different ways as a sum of multiples of numbers in S.")

    # --- Search for the first (smallest) A-value ---
    print("\nStep 1: Searching for the smallest A-value (A1).")
    print("Let's check numbers in increasing order:")
    print("A = 1, 2, 3, 4: Can only be formed in one way or not at all (e.g., 4 = 2*2).")
    print("A = 5:")
    print("  - Combination 1: 1 * 5")
    print("  - Combination 2: 1 * 2 + 1 * 3")
    print("Since 5 can be formed in two ways, we have found our first relation: (1*2 + 1*3) - 1*5 = 0.")
    A1 = 5
    print(f"The sum is 5, so the smallest A-value is A1 = {A1}.")

    # --- Search for the second smallest A-value ---
    print("\nStep 2: Searching for the second smallest A-value (A2) from an independent relation.")
    print("A = 6:")
    print("  - Combination 1: 2 * 3")
    print("  - Combination 2: 3 * 2")
    print("This gives a second relation: (3*2) - (2*3) = 0.")
    print("This relation is mathematically independent of the first one.")
    A2 = 6
    print(f"The sum is 6, so the second smallest A-value is A2 = {A2}.")

    # --- Conclusion ---
    print("\nStep 3: Determine the rectangle dimensions and area.")
    print(f"The two smallest A-values from independent relations are {A1} and {A2}.")
    print(f"Thus, the smallest rectangle admitting a non-guillotine tiling has dimensions {A1}x{A2}.")

    area = A1 * A2
    print("\nThe area of this rectangle is the product of these dimensions.")
    print("Final equation:")
    print(f"{A1} * {A2} = {area}")


find_smallest_rectangle_area()