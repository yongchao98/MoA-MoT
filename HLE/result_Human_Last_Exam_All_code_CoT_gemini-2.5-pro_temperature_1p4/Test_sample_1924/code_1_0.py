def solve_geometry_problem(MA, MB):
    """
    Calculates the value of MG - MH based on the theorem that MG - MH = MA - MB.

    Args:
        MA (float): The length of the segment MA.
        MB (float): The length of the segment MB.
    """
    # According to the theorem for this geometric configuration, MG - MH = MA - MB.
    MG_minus_MH = MA - MB

    # Print the final equation with the given numbers
    print(f"Given MA = {MA} and MB = {MB}")
    print(f"The expression for MG - MH is MA - MB.")
    print(f"Therefore, MG - MH = {MA} - {MB} = {MG_minus_MH}")

# Example usage with arbitrary values for MA and MB.
# Let's assume M is a point on the chord AB, with MA=7 and MB=3.
MA = 7.0
MB = 3.0
solve_geometry_problem(MA, MB)
