import collections

def solve():
    """
    This function encapsulates the reasoning for classifying the rings and prints the final result.
    """
    # Based on mathematical analysis, the rings are partitioned into the following isomorphism classes:
    
    # Class 1: {A, B}
    # Rings A and B are coordinate rings of elliptic curves defined by y^2 = x^3 + x^2 - 3x + 1
    # and y^2 = x^3 + 2x^2 - 2x + 3, respectively. A change of variables x' = x + 5
    # transforms the first equation into the second (modulo 7), showing they are isomorphic.
    # They are infinite integral domains, but not local rings.
    
    # Class 2: {C, E, L}
    # Ring C is F_7[x]/(5x^2+x+1). The polynomial factors as 5(x-1)(x-3).
    # Ring E is F_7[x]/(3x^2+x+6). The polynomial factors as 3x(x-2).
    # By the Chinese Remainder Theorem, both C and E are isomorphic to F_7 x F_7.
    # Ring L is defined as F_7 x F_7. Thus, C, E, and L are isomorphic.
    
    # Class 3: {D, H}
    # Ring D is defined by an ideal whose Gröbner basis over F_7 is {1}. This means the ideal
    # is the entire ring F_7[x,y], and the quotient ring D is the zero ring {0}.
    # Ring H is a quotient of F_7[[x]]. The generator of the ideal is a unit because its
    # constant term is 4/4 = 1. Thus, the ideal is F_7[[x]] itself, and H is also the zero ring.
    
    # Class 4: {F, G}
    # Ring F is F_7[x]/(x^2).
    # Ring G is F_7[x]/(x^2+3x+4). The polynomial is a perfect square: (x-2)^2.
    # The map t -> x-2 gives an isomorphism from F_7[t]/(t^2) to G.
    # These rings are local, have size 49, and contain non-zero nilpotent elements.
    
    # Class 5: {I}
    # Ring I is the coordinate ring of the elliptic curve y^2 = x^3+3x^2+3x+2 = (x+1)^3+1.
    # This curve is not isomorphic to the curve corresponding to A and B over F_7.
    # It is an infinite, non-local integral domain, distinct from all other rings.
    
    # Class 6: {J}
    # Ring J is the local ring of the affine line at the point (x+1), which is a Discrete Valuation Ring (DVR).
    # As a DVR, it is a local ring. It is an infinite integral domain but not a finitely-generated
    # F_7-algebra, which distinguishes it from A, B, and I.
    
    # Class 7: {K}
    # Ring K is the finite field F_49. It is the only ring in the list that is a field of size 49.
    
    # Combine the letters for each class, sort them alphabetically.
    classes = {
        'A': 'AB', 'B': 'AB',
        'C': 'CEL', 'E': 'CEL', 'L': 'CEL',
        'D': 'DH', 'H': 'DH',
        'F': 'FG', 'G': 'FG',
        'I': 'I',
        'J': 'J',
        'K': 'K'
    }
    
    # Use an ordered dictionary to get unique class strings, preserving insertion order.
    # Sorting keys of the dictionary to ensure alphabetical order of groups.
    unique_classes = collections.OrderedDict()
    for key in sorted(classes.keys()):
        unique_classes[classes[key]] = True
        
    result = "[" + ", ".join(unique_classes.keys()) + "]"
    
    print(result)

solve()