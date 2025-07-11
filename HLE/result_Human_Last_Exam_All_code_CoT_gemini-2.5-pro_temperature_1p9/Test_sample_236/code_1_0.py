import textwrap

def solve_homology_cobordism_question():
    """
    Calculates the number of homology cobordism group elements representable
    by integral surgery on knots with at most four crossings.

    The logic relies on established results from low-dimensional topology.
    """
    print(textwrap.dedent("""\
    Analyzing integral surgeries on knots with <= 4 crossings to find the number of
    representable elements in the homology cobordism group.

    The knots with at most 4 crossings are:
    - 0_1 (the unknot)
    - 3_1 (the trefoil knot)
    - 4_1 (the figure-eight knot)

    For the result of an integral p/q surgery on a knot K to be a homology sphere,
    we must have |p|=1. So we only consider +1 and -1 surgery.
    """))

    # We use symbolic names for the homology cobordism classes.
    # 'id': The identity element (S^3).
    # 'A': The class from +1 surgery on the right-handed trefoil.
    # 'B': The class of the Poincaré homology sphere.
    # '-A', '-B': The inverses of 'A' and 'B'.
    
    all_elements = set()
    
    # 1. Unknot (0_1) - Amphichiral
    # +/-1 surgery on the unknot yields S^3, the identity element.
    num_before = len(all_elements)
    all_elements.add('id')
    num_from_unknot = len(all_elements) - num_before
    print(textwrap.dedent(f"""\
    -----------------------------------------------------
    Knot: Unknot (0_1)
    - Surgeries: +1, -1
    - Results: Both surgeries yield the 3-sphere S^3, which is the identity element ('id').
    - New elements contributed: {num_from_unknot}
    """))
    
    # 2. Trefoil (3_1) - Chiral
    # We must consider the right-handed (RHT) and left-handed (LHT) trefoils.
    # Let +1 surgery on RHT be class 'A' (this is Sigma(2,3,7)).
    # Let -1 surgery on RHT be class 'B' (this is the Poincaré sphere).
    # Surgery on the mirror image LHT gives inverse classes:
    # M(LHT, p) is the inverse of M(RHT, -p).
    # So, +1 on LHT gives -B, and -1 on LHT gives -A.
    # Topological Fact: The Poincaré sphere has order 2, so B = -B.
    num_before = len(all_elements)
    trefoil_elements = {'A', 'B', '-A'} # B and -B are the same element
    all_elements.update(trefoil_elements)
    num_from_trefoil = len(all_elements) - num_before
    print(textwrap.dedent(f"""\
    -----------------------------------------------------
    Knot: Trefoil (3_1) and its mirror
    - R-Trefoil surgeries: +1 gives class 'A', -1 gives class 'B' (Poincaré sphere).
    - L-Trefoil surgeries: +1 gives '-B', -1 gives '-A'.
    - Applying identity B = -B (Poincaré sphere has order 2), the new distinct classes are 'A', '-A', and 'B'.
    - New elements contributed: {num_from_trefoil}
    """))
    
    # 3. Figure-eight (4_1) - Amphichiral
    # It's known that -1 surgery on 4_1 gives the Poincaré sphere (class 'B').
    # Since 4_1 is amphichiral, +1 surgery gives the inverse element, -B.
    # Since B = -B, both surgeries give the same element 'B'.
    num_before = len(all_elements)
    fig8_elements = {'B'}
    all_elements.update(fig8_elements)
    num_from_fig8 = len(all_elements) - num_before
    print(textwrap.dedent(f"""\
    -----------------------------------------------------
    Knot: Figure-eight (4_1)
    - Surgeries: +1, -1
    - Results: -1 surgery gives the Poincaré sphere ('B'). +1 surgery gives '-B'.
    - Since B = -B, both surgeries produce the same element 'B', which was already found.
    - New elements contributed: {num_from_fig8}
    """))
    
    total_count = len(all_elements)
    
    print("-----------------------------------------------------")
    print("\nTotal Count Calculation:")
    print("The final count is the sum of new, unique elements from each knot type.")
    
    final_equation = f"{num_from_unknot} (from unknot) + {num_from_trefoil} (from trefoil) + {num_from_fig8} (from figure-eight) = {total_count}"
    
    print("\n" + final_equation)

if __name__ == '__main__':
    solve_homology_cobordism_question()
<<<4>>>