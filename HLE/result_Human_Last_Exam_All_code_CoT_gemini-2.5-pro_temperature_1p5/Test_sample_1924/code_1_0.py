import sympy

def solve_geometry_problem():
    """
    Solves the geometry problem symbolically.

    The problem describes a geometric configuration involving a circle O, a chord AB,
    and a point M on AB. Two other chords, CD and EF, are used to define two new circles:
    the circumcircle of triangle EMD and the circumcircle of triangle CMF. These
    circles intersect the line AB at G and H, respectively (apart from M).
    The goal is to express the difference MG - MH in terms of MA and MB.

    This problem is a known result in geometry, sometimes called Stevanovic's Theorem.
    The theorem states that for the signed lengths along the line AB, the following
    relationship holds:
    vector(MG) - vector(MH) = vector(MB) - vector(MA)

    The quantities MA and MB in the problem are lengths, which are positive.
    The expression MG - MH refers to the difference of signed lengths.

    Let's set up a coordinate system on the line AB, with M at the origin.
    Let the direction from A to B be positive.
    If M is between A and B:
      - The coordinate of A is -MA.
      - The coordinate of B is +MB.
    The signed length of a segment from M to a point X, denoted as MX, is simply
    the coordinate of X.
    So, the signed length MA is -MA (the length) and the signed length MB is +MB.
    
    The theorem gives:
    MG - MH = (signed MB) - (signed MA) = (+MB) - (-MA) = MA + MB.
    This result (MA+MB = AB) is one interpretation.

    However, another common statement of the theorem is MA-MB = MH-MG. If we use this
    statement with signed lengths (let's denote signed XY as xy):
    ma - mb = mh - mg
    Rearranging this gives:
    mg - mh = mb - ma
    
    This does not depend on M's position relative to A and B. It is a more robust result.
    Let's express the result in terms of the lengths MA and MB. The expression 'mb - ma'
    is symbolic and its numerical value depends on the relative positions of A, M, B.
    Without knowing the numerical values or the exact positions, the most precise answer
    is the symbolic expression itself.
    
    The code below will represent MA and MB as symbolic variables and print the
    resulting expression for MG - MH.
    """
    
    MA, MB = sympy.symbols('MA MB')
    
    # According to the theorem, the difference of signed lengths MG - MH is equal to MB - MA.
    result_expression = MB - MA
    
    print("The value of MG - MH is expressed in terms of the lengths MA and MB.")
    # We print the result symbolically. The numbers in the final equation are 1 and -1,
    # the coefficients of MB and MA.
    print(f"MG - MH = (1)*MB + (-1)*MA")
    print(f"So, MG - MH = {result_expression}")

solve_geometry_problem()