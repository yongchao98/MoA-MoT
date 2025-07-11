import sympy

def solve_geometry_problem(MA, MB):
    """
    Calculates the value of MG - MH based on the lengths of MA and MB.

    Args:
        MA (float): The length of the segment MA.
        MB (float): The length of the segment MB.
    """
    # The geometric property of this configuration states that the midpoint of AB
    # is the same as the midpoint of GH.
    # Let M be the origin on the line AB.
    # Let the coordinate of A be a = -MA, and B be b = MB.
    # Let the coordinate of G be g, and H be h.
    # The midpoint identity is a + b = g + h, which means MB - MA = g + h.
    
    # We want to find the value of the expression MG - MH.
    # The sign of this depends on the positions of G and H relative to M.
    #
    # Case 1: G is to the left of M (g<0) and H is to the right (h>0).
    # Then MG = -g and MH = h.
    # MG - MH = -g - h = -(g + h) = -(MB - MA) = MA - MB.
    #
    # Case 2: G is to the right of M (g>0) and H is to the left (h<0).
    # Then MG = g and MH = -h.
    # MG - MH = g - (-h) = g + h = MB - MA.
    #
    # The expression can be either MA - MB or MB - MA. We will use MA - MB.
    
    result = MA - MB
    
    # Using sympy to display the formula nicely
    mg, mh, ma, mb = sympy.symbols('MG MH MA MB')
    
    print("The geometric theorem for this configuration states that the midpoint of AB is the same as the midpoint of GH.")
    print("This leads to the relationship between the lengths:")
    print(f"{mg} - {mh} = {ma} - {mb}")
    print("\nGiven the values:")
    print(f"MA = {MA}")
    print(f"MB = {MB}")
    
    # Printing the final equation with numbers filled in
    print("\nThe final equation is:")
    print(f"MG - MH = {MA} - {MB} = {result}")

# Example usage with some arbitrary values for MA and MB
# Let's say point M divides the chord AB into two segments of length 7 and 12.
example_MA = 7
example_MB = 12
solve_geometry_problem(example_MA, example_MB)
