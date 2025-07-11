def solve_geometry_ratio():
    """
    This function derives the expression for the ratio BM/MI in terms of triangle side lengths a, b, and c.
    The derivation is explained step-by-step using print statements.
    """
    print("Let a, b, and c be the lengths of sides BC, CA, and AB, respectively.")
    print("Let the angles of the triangle be α, β, and γ at vertices A, B, and C.")
    print("-" * 50)

    print("Step 1: Apply Ptolemy's Theorem.")
    print("The points A, M, C, and B lie on the circumcircle, forming a cyclic quadrilateral AMCB.")
    print("Ptolemy's Theorem states that for a cyclic quadrilateral, the sum of the products of opposite sides equals the product of the diagonals.")
    print("Applying this to AMCB, we get:")
    print("AM * CB + CM * AB = AC * BM")
    print("Substituting the side lengths a, b, c:")
    print("AM * a + CM * c = b * BM")
    print("-" * 50)

    print("Step 2: Use the property that M is on the angle bisector of ∠B.")
    print("Since M is on the circumcircle, the chords AM and CM subtend the angles ∠ABM and ∠CBM at the circumference.")
    print("As BM is the angle bisector of ∠B, we have ∠ABM = ∠CBM = β/2.")
    print("Chords that subtend equal angles at the circumference are equal in length. Therefore, AM = CM.")
    print("Substituting AM = CM into the equation from Step 1:")
    print("AM * a + AM * c = b * BM")
    print("AM * (a + c) = b * BM")
    print("Rearranging this gives us a ratio for BM / AM:")
    print("BM / AM = (a + c) / b")
    print("-" * 50)

    print("Step 3: Prove that AM = MI.")
    print("This is a well-known property. Let's prove it by considering the angles in triangle AIM.")
    print("Angle ∠MAI = ∠MAC + ∠CAI.")
    print(" - ∠MAC subtends the same arc (MC) as ∠MBC, so ∠MAC = ∠MBC = β/2.")
    print(" - AI is the angle bisector of ∠A, so ∠CAI = α/2.")
    print(" - Therefore, ∠MAI = α/2 + β/2.")
    print("Angle ∠MIA is an exterior angle to triangle ABI.")
    print(" - Therefore, ∠MIA = ∠IAB + ∠IBA = α/2 + β/2.")
    print("Since ∠MAI = ∠MIA, triangle AIM is an isosceles triangle with the sides opposite to these angles being equal.")
    print("So, AM = MI.")
    print("-" * 50)

    print("Step 4: Combine the results to find the final ratio.")
    print("From Step 2, we have the ratio: BM / AM = (a + c) / b")
    print("From Step 3, we have the equality: AM = MI")
    print("Substitute MI for AM in the ratio:")
    print("BM / MI = (a + c) / b")
    print("-" * 50)

    print("The final expression for the ratio is:")
    print("BM   a + c")
    print("-- = -----")
    print("MI     b")

solve_geometry_ratio()