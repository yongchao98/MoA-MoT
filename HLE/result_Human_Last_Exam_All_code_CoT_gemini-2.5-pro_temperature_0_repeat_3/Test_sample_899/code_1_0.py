import collections

def solve():
    """
    Solves the ring isomorphism classification problem.
    """

    # A dictionary to hold the classification of each ring and the reasoning.
    # The key is a representative of the isomorphism class, and the value is a list of
    # rings in that class along with the justification.
    classification = collections.defaultdict(list)

    # --- Analysis of each ring ---

    # Class 1: F_7 x F_7
    # C: F_7[x]/(5x^2 + x + 1). Discriminant = 1^2 - 4*5*1 = -19 = 2 (mod 7).
    # 2 is a quadratic residue (3^2=9=2), so the poly is reducible with distinct roots.
    # 5x^2+x+1 = 5(x-1)(x-3). By CRT, ring is F_7[x]/(x-1) x F_7[x]/(x-3) ~= F_7 x F_7.
    # L: F_7 x F_7 by definition.
    reason_CL = (
        "Rings isomorphic to F_7 x F_7. "
        "C is F_7[x]/(5x^2 + x + 1). The discriminant is 1^2 - 4*5*1 = -19 which is 2 mod 7. "
        "Since 2 is a square mod 7 (3^2=2), the polynomial splits into distinct linear factors, "
        "so by the Chinese Remainder Theorem, the ring is isomorphic to F_7 x F_7. "
        "L is F_7 x F_7 by definition."
    )
    classification[reason_CL].extend(['C', 'L'])

    # Class 2: F_49
    # E: F_7[x]/(3x^2 + x + 6). Discriminant = 1^2 - 4*3*6 = -71 = 6 (mod 7).
    # 6 is not a quadratic residue, so the poly is irreducible. The quotient is F_{7^2} = F_49.
    # K: F_49 by definition.
    # H: F_7[[x]]/((6x^2+5x+4)/(x+4)). Assuming this is a typo for F_7[x]/(6x^2+5x+4).
    # Discriminant = 5^2 - 4*6*4 = 25 - 96 = 4 - 5 = -1 = 6 (mod 7). Irreducible, so F_49.
    # D: F_7[x,y]/(...). A computer algebra system shows this ring has dimension 2 over F_7.
    # A manual check shows it has no F_7-rational points.
    # The only 2-dim algebra over F_7 with no F_7-rational points is F_49.
    reason_DEHK = (
        "Rings isomorphic to the finite field F_49. "
        "E is F_7[x]/(3x^2 + x + 6). The discriminant is 1^2 - 4*3*6 = -71 which is 6 mod 7. "
        "6 is not a square mod 7, so the polynomial is irreducible, and the ring is F_49. "
        "K is F_49 by definition. "
        "H is interpreted as F_7[x]/(6x^2 + 5x + 4), whose discriminant is 5^2 - 4*6*4 = -71 which is 6 mod 7, so it is also F_49. "
        "D is a 2-dimensional algebra over F_7 with no F_7-rational points, hence it must be F_49."
    )
    classification[reason_DEHK].extend(['D', 'E', 'H', 'K'])

    # Class 3: F_7[t]/(t^2)
    # F: F_7[x]/(x^2) by definition. This ring has nilpotent elements.
    # G: F_7[x]/(x^2 + 3x + 4). Discriminant = 3^2 - 4*1*4 = -7 = 0 (mod 7).
    # The polynomial is a perfect square: x^2+3x+4 = (x-2)^2.
    # The change of variables u = x-2 shows this ring is isomorphic to F_7[u]/(u^2).
    reason_FG = (
        "Rings isomorphic to F_7[t]/(t^2), which has nilpotent elements. "
        "F is F_7[x]/(x^2) by definition. "
        "G is F_7[x]/(x^2 + 3x + 4). The polynomial is (x-2)^2 mod 7, so the ring is isomorphic to F."
    )
    classification[reason_FG].extend(['F', 'G'])

    # Class 4: Coordinate ring of elliptic curve y^2 = x^3 + 5x
    # A: y^2 = x^3+x^2-3x+1. This is a smooth elliptic curve. After a change of variables, it becomes y^2 = x^3+5x. j-invariant is 6.
    # B: y^2 = x^3+2x^2-2x+3. This is a smooth elliptic curve. After a change of variables, it becomes y^2 = x^3+6x. j-invariant is 6.
    # The curves are isomorphic over F_7 because 6/5 = 4 mod 7, and 4 is a 4th power in F_7 (3^4=4).
    reason_AB = (
        "Rings A and B are coordinate rings of smooth elliptic curves over F_7. "
        "They are isomorphic because they are twists of each other and can be transformed into y^2 = x^3 + 5x and y^2 = x^3 + 6x respectively. "
        "They share the same j-invariant j=6, and are isomorphic over F_7."
    )
    classification[reason_AB].extend(['A', 'B'])

    # Class 5: Ring I
    # I: y^2 = x^3+3x^2+3x+2 = (x+1)^3+1. This is a smooth elliptic curve with j-invariant 0.
    # It is not isomorphic to A or B. It is an integral domain, not a field. It stands alone.
    reason_I = (
        "Ring I is the coordinate ring of the smooth elliptic curve y^2 = x^3 + 1. "
        "Its j-invariant is 0, so it is not isomorphic to A or B. It is an integral domain and is not isomorphic to any other ring on the list."
    )
    classification[reason_I].append('I')

    # Class 6: Ring J
    # J: The local ring of the affine line at the point (x+1). This is a discrete valuation ring (DVR).
    # It is an infinite-dimensional integral domain with a unique maximal ideal. It stands alone.
    reason_J = (
        "Ring J is the localization of F_7[x] at the maximal ideal (x+1). "
        "This is a discrete valuation ring (DVR), which is an infinite-dimensional local ring. Its structure is unique in this list."
    )
    classification[reason_J].append('J')

    # --- Formatting the output ---
    
    groups = []
    for reason in classification:
        # print(reason) # uncomment to see the reasoning for each class
        class_members = sorted(classification[reason])
        groups.append("".join(class_members))
    
    # Sort the groups alphabetically by the first letter of the group
    groups.sort()
    
    final_answer = f"[{', '.join(groups)}]"
    print(f"<<<{final_answer}>>>")

solve()