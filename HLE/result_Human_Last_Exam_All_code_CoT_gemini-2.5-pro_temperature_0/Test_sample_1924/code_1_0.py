def solve_geometry_problem(MA, MB):
    """
    Calculates the value of MG - MH based on the lengths of MA and MB.

    In the given geometric configuration, there is a remarkable relationship between
    the segments on the line AB. The relationship is:
    MG - MH = MB - MA

    This expression holds true for the lengths of the segments, where the
    differences can be positive or negative.

    Args:
        MA (float): The length of the segment MA.
        MB (float): The length of the segment MB.
    """
    # The problem is to express MG - MH in terms of MA and MB.
    # The direct relationship is MG - MH = MB - MA.
    result = MB - MA

    print(f"Given lengths are:")
    print(f"MA = {MA}")
    print(f"MB = {MB}")
    print("\nThe relationship is: MG - MH = MB - MA")
    print("\nCalculating the value of the expression:")
    # We output each number in the final equation.
    # Since we don't know MG and MH individually, we calculate the final result.
    print(f"MG - MH = {MB} - {MA} = {result}")

# Example usage with some arbitrary values for MA and MB.
# Let's assume M is a point on a chord AB of length 15,
# such that MA = 10 and MB = 5.
MA_example = 10.0
MB_example = 5.0
solve_geometry_problem(MA_example, MB_example)
