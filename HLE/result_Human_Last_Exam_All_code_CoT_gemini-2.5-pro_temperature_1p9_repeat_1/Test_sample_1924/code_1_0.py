import math

def solve_geometry_problem(MA, MB):
    """
    Calculates the value of MG - MH based on the lengths of MA and MB.

    In the given geometric configuration, the following relationship holds:
    MG - MH = MA - MB

    Args:
        MA (float): The length of the segment MA.
        MB (float): The length of the segment MB.

    Returns:
        float: The value of MG - MH.
    """
    # According to the geometric theorem related to this construction,
    # the expression MG - MH is equal to MA - MB.
    result = MA - MB
    
    # We print the final equation with the substituted values.
    # Note: MA, MB, MG, MH can be considered as directed segments (vectors) on the line AB.
    # Here, MA and MB are given as lengths (magnitudes). The expression MA - MB gives the
    # value of MG - MH, which is also a directed length.
    print(f"Given MA = {MA}")
    print(f"Given MB = {MB}")
    print(f"The relationship is MG - MH = MA - MB")
    print(f"Therefore, MG - MH = {MA} - {MB} = {result}")

# You can change these example values for MA and MB.
# For example, let's say M divides the chord AB such that MA = 7 and MB = 3.
MA_val = 7.0
MB_val = 3.0

solve_geometry_problem(MA_val, MB_val)
