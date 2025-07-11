def calculate_periods():
    """
    Calculates the period for each case based on the number of vertices (N)
    and the order of rotational symmetry (g) of the polygon formed by the points.
    The period is calculated as N/g.
    """

    # Case 1: 13, 31, 23
    N1 = 3  # Equilateral triangle
    g1 = 3  # C3 symmetry
    period1 = N1 // g1
    print(f"Case 1: Points {13, 31, 23} form a shape with {N1} vertices.")
    print(f"The shape is an equilateral triangle with 3-fold rotational symmetry (g={g1}).")
    print(f"Period = {N1} / {g1} = {period1}")
    print("-" * 20)

    # Case 2: 10, 4, 23, 31
    N2 = 4  # Rectangle
    g2 = 2  # C2 symmetry
    period2 = N2 // g2
    print(f"Case 2: Points {10, 4, 23, 31} form a shape with {N2} vertices.")
    print(f"The shape is a rectangle with 2-fold rotational symmetry (g={g2}).")
    print(f"Period = {N2} / {g2} = {period2}")
    print("-" * 20)

    # Case 3: 5, 15, 17, 19, 21, 7
    N3 = 6  # Asymmetric polygon
    g3 = 1  # C1 symmetry (no symmetry)
    period3 = N3 // g3
    print(f"Case 3: Points {5, 15, 17, 19, 21, 7} form a shape with {N3} vertices.")
    print(f"The shape is assumed to be asymmetric, with 1-fold rotational symmetry (g={g3}).")
    print(f"Period = {N3} / {g3} = {period3}")
    print("-" * 20)

    # Case 4: 4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13
    N4 = 13 # Asymmetric polygon
    g4 = 1   # C1 symmetry (no symmetry)
    period4 = N4 // g4
    print(f"Case 4: Points from the long list form a shape with {N4} vertices.")
    print(f"The shape is assumed to be asymmetric, with 1-fold rotational symmetry (g={g4}).")
    print(f"Period = {N4} / {g4} = {period4}")
    print("-" * 20)
    
    # Final answer string
    final_answer = f"{period1},{period2},{period3},{period4}"
    print(f"The four periods separated by ',' are: {final_answer}")

calculate_periods()