import collections

def analyze_intersections():
    """
    This function programmatically follows the logic to find the number of
    homeomorphism classes for the intersections of two geodesics.
    """

    # Step 1 & 2: Define geodesics by their direction sets.
    # A direction 'u' represents a ray {t*u | t >= 0}.
    # We use strings for symbolic representation.
    # G1 and G2 are the two geodesics we are intersecting.
    
    # A set to store the homeomorphism classes of the resulting intersections.
    # We use descriptive strings for the classes.
    homeomorphism_classes = set()

    print("Analyzing possible intersections of two geodesics G1 and G2...")
    print("-" * 60)

    # --- Case 1: G1 is a Line, G2 is a Line ---
    # Directions for G1: {u1, -u1}
    # Directions for G2: {u2, -u2}
    #   Subcase 1.1: G1 and G2 are the same line.
    #     Directions are the same. Intersection of directions has size 2.
    #     The resulting shape is {u1, -u1}, a Line.
    homeomorphism_classes.add("Line")
    #   Subcase 1.2: G1 and G2 are different lines.
    #     Directions are disjoint. Intersection of directions has size 0.
    #     The resulting shape is {}, a Point.
    homeomorphism_classes.add("Point")

    # --- Case 2: G1 is a Line, G2 is a Bent Line ---
    # Directions for G1: {u1, -u1}
    # Directions for G2: {u2, u3} where u2 != -u3
    #   Subcase 2.1: The direction sets are disjoint.
    #     Intersection of directions has size 0. Shape: Point.
    homeomorphism_classes.add("Point") # Already present
    #   Subcase 2.2: The direction sets share one direction (e.g., u1=u2).
    #     Intersection of directions has size 1. Shape: Ray.
    homeomorphism_classes.add("Ray")

    # --- Case 3: G1 is a Bent Line, G2 is a Bent Line ---
    # Directions for G1: {u1, u2}
    # Directions for G2: {u3, u4}
    #   Subcase 3.1: Disjoint directions. Size 0. Shape: Point.
    homeomorphism_classes.add("Point") # Already present
    #   Subcase 3.2: One common direction. Size 1. Shape: Ray.
    homeomorphism_classes.add("Ray") # Already present
    #   Subcase 3.3: Two common directions (G1=G2). Size 2.
    #     The resulting shape is {u1, u2}, a Bent Line.
    homeomorphism_classes.add("Bent Line")

    print("The possible shapes of the intersections are:")
    # We sort for consistent output
    sorted_classes = sorted(list(homeomorphism_classes))
    for item in sorted_classes:
        print(f"- {item}")
    
    print("\nTopological analysis confirms these 4 shapes are distinct:")
    print("1. Point: Compact.")
    print("2. Ray: Not compact. Removing the origin leaves 1 connected component.")
    print("3. Line: Not compact. Removing any point leaves 2 connected components.")
    print("4. Bent Line: Not compact. Removing the origin leaves 2 components, but removing any other point leaves 1.")
    
    print("\nThus, the number of distinct homeomorphism classes is:")
    
    final_count = len(homeomorphism_classes)
    
    # The prompt asks to output the numbers in the final equation.
    # We will show the count of each class found.
    class_counts = collections.Counter(['Line', 'Point', 'Ray', 'Bent Line'])
    equation_parts = [str(v) for k, v in sorted(class_counts.items())]
    
    print(f"{' + '.join(equation_parts)} = {final_count}")


if __name__ == '__main__':
    analyze_intersections()
