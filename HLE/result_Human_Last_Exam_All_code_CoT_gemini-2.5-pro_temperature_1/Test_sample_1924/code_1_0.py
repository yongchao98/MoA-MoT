def solve_geometry_problem(ma, mb):
    """
    This function calculates the value of MG - MH based on the lengths of MA and MB.

    As derived from the geometric properties of the setup, the relationship is:
    MG - MH = MA - MB

    Args:
        ma (float): The length of the segment MA.
        mb (float): The length of the segment MB.
    """
    
    # The value of MG - MH is simply MA - MB.
    result = ma - mb
    
    print("The relationship is: MG - MH = MA - MB")
    print(f"Given MA = {ma} and MB = {mb}:")
    
    # Print the equation with the given numbers
    print(f"MG - MH = {ma} - {mb}")
    
    # Print the final result
    print(f"MG - MH = {result}")

# You can change these example values for MA and MB
MA_length = 15
MB_length = 7

solve_geometry_problem(MA_length, MB_length)