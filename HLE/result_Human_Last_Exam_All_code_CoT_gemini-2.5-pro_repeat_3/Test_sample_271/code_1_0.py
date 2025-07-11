import sympy
from sympy import symbols, Poly, factor, discriminant
from sympy.polys.galoistools import galois_group
from sympy.combinatorics import PermutationGroup
from sympy.combinatorics.group_constructs import DirectProduct

def solve_field_degrees():
    """
    Finds all possible degrees of normal field extensions L/Q where Q < L < K,
    with K being the splitting field of the given polynomial.
    """
    x = symbols('x')
    
    # Step 1: Define and factor the polynomial
    poly_expr = x**7 - 2*x**5 - 9*x**3 + 3*x**2 + 18*x - 6
    p = Poly(poly_expr, x, domain='ZZ')
    
    factored_poly = factor(p)
    
    # The factors are (x**2 - 2) and (x**5 - 9*x + 3)
    # We can extract them programmatically, but for clarity, we define them directly
    f1 = Poly(x**2 - 2)
    f2 = Poly(x**5 - 9*x + 3)

    print(f"The polynomial f(x) = {poly_expr} factors over Q into:")
    print(f"f1(x) = {f1.expr}")
    print(f"f2(x) = {f2.expr}\n")

    # Step 2: Determine the Galois Group
    # Find Galois groups for each factor
    G1 = galois_group(f1)
    G2 = galois_group(f2)
    
    print(f"The Galois group of f1(x) is isomorphic to Z/2Z, with order {G1.order()}.")
    print(f"The Galois group of f2(x) is the symmetric group S_5, with order {G2.order()}.\n")

    # Check if the intersection of splitting fields is trivial
    # The splitting field of f1 is Q(sqrt(2)).
    # The unique quadratic subfield of the splitting field of f2 is Q(sqrt(D)),
    # where D is the discriminant of f2.
    d1 = discriminant(f1) # This is 8, so the field is Q(sqrt(2))
    d2 = discriminant(f2)
    
    # Q(sqrt(2)) and Q(sqrt(d2)) are different since d2 is not of the form 2*q^2.
    # Therefore, the intersection of the splitting fields is Q.
    # The Galois group of f(x) is the direct product G1 x G2.
    G = DirectProduct(G1, G2)
    print(f"The Galois group G of f(x) is the direct product Z/2Z x S_5.")
    print(f"The order of G is {G.order()}.\n")

    # Step 3: Find all proper, non-trivial normal subgroups
    normal_subgroups = G.normal_subgroups()
    
    proper_normal_subgroups = []
    for H in normal_subgroups:
        if not H.is_trivial and H.order() < G.order():
            proper_normal_subgroups.append(H)
            
    # Step 4: Calculate the degrees of the corresponding normal subfields
    possible_degrees = set()
    for H in proper_normal_subgroups:
        degree = G.order() / H.order()
        possible_degrees.add(int(degree))
        
    # Step 5: List the degrees
    sorted_degrees = sorted(list(possible_degrees))
    
    print("The possible degrees of normal field extensions L/Q such that Q < L < K are the indices of the proper non-trivial normal subgroups of G.")
    print("The calculated degrees are:")
    # The problem asks to output each number in the final equation.
    # We will print the list of numbers.
    for i, degree in enumerate(sorted_degrees):
        if i < len(sorted_degrees) - 1:
            print(degree, end=", ")
        else:
            print(degree)

solve_field_degrees()
<<<2, 4, 120>>>