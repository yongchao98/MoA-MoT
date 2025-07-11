def solve_geometry_ratio():
    """
    This function solves the geometry problem by explaining the derivation steps
    and printing the final expression for the ratio BM / MI.
    """
    
    # Define the side lengths as strings for display purposes
    a = 'a'
    b = 'b'
    c = 'c'
    
    # Print the explanation of the derivation
    print("This program solves for the ratio BM/MI in terms of the side lengths a, b, and c of triangle ABC.")
    print("The derivation follows these steps:\n")

    print("Step 1: Express the ratio BM/MI in a simpler form.")
    print("Since B, I, and M are collinear, BM = BI + IM.")
    print("Therefore, BM / MI = (BI + IM) / IM = (BI / IM) + 1.\n")

    print("Step 2: Use the Trillium Theorem (or Incenter-Excenter Lemma).")
    print("This theorem states that the point M where the angle bisector from B intersects the circumcircle is equidistant from A, C, and the incenter I.")
    print("Thus, we have the crucial identity: IM = MA.\n")
    
    print("Step 3: Substitute and simplify the problem.")
    print("Replacing IM with MA, the ratio becomes (BI / MA) + 1.")
    print("The problem is now to find the ratio BI / MA.\n")
    
    print("Step 4: Express BI and MA in terms of angles and the circumradius R.")
    print("Using the Law of Sines and properties of the triangle, we find:")
    print("  BI = 2R * sin(A/2) * sin(C/2) / sin(B/2)")
    print("  MA = 2R * sin(B/2)\n")

    print("Step 5: Calculate the ratio BI / MA and convert it to side lengths.")
    print("The ratio BI / MA simplifies to 2 * sin(A/2) * sin(C/2) / sin(B/2).")
    print("Using the half-angle formulas to express this in terms of side lengths a, b, and c gives:")
    print(f"  BI / MA = ({a} + {c} - {b}) / {b}\n")
    
    print("Step 6: Compute the final result.")
    print("Substitute the expression for BI / MA back into our ratio from Step 3:")
    print(f"  BM / MI = (BI / MA) + 1 = (({a} + {c} - {b}) / {b}) + 1")
    print(f"  BM / MI = ({a} + {c} - {b} + {b}) / {b}\n")

    print("-------------------------------------------")
    print("Final Expression:")
    # The instruction "output each number in the final equation!" is interpreted
    # as printing the variables in the final formula.
    print(f"The ratio BM / MI is equal to the expression: ({a} + {c}) / {b}")
    print("-------------------------------------------")

# Execute the function to display the solution
solve_geometry_ratio()
