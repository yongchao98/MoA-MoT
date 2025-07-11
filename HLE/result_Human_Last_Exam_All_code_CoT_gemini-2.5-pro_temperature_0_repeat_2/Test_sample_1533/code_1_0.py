def solve_geometry_ratio():
    """
    This function prints the derived expression for the ratio BM/MI.
    
    Let a, b, c be the lengths of the sides of triangle ABC opposite to
    vertices A, B, and C, respectively.
    a = BC
    b = AC
    c = AB
    
    The ratio BM/MI is derived using the Incenter-Excenter Lemma and Ptolemy's Theorem.
    """
    
    # Define the side lengths symbolically for printing
    a = 'a'
    b = 'b'
    c = 'c'
    
    # The final derived equation is BM/MI = (a + c) / b
    numerator = f"({a} + {c})"
    denominator = f"{b}"
    
    print("The expression for the ratio BM/MI in terms of the side lengths a, b, and c is:")
    print(f"  BM   {numerator}")
    print(f"---- = {'''-''' * (len(numerator) + 2)}")
    print(f"  MI    {denominator}")

solve_geometry_ratio()