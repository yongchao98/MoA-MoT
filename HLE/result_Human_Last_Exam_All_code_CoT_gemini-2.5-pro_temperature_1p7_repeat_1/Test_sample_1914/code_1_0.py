# This script calculates the number of categories with 2 objects and 4 morphisms, up to isomorphism.
# The calculation is based on a case-by-case analysis of how the 4 morphisms can be distributed.

# A category with 2 objects (A, B) and 4 morphisms requires:
# - An identity morphism for each object: id_A in Hom(A,A) and id_B in Hom(B,B).
# - The total number of morphisms |Hom(A,A)| + |Hom(B,B)| + |Hom(A,B)| + |Hom(B,A)| = 4.

# We analyze the problem by considering the non-isomorphic distributions of morphisms.
# The number of categories is the sum of counts from each distinct structural pattern.

# Pattern 1: Morphism distribution is (3, 1, 0, 0) or (1, 3, 0, 0).
# In this case, one object has 3 endomorphisms (forming a monoid of order 3),
# and the other has only the identity.
# The number of non-isomorphic monoids of order 3 is a known result from algebra, which is 5.
# These 5 monoids define 5 distinct categories.
count_pattern1 = 5

# Pattern 2: Distribution is (1, 1, 2, 0) or (1, 1, 0, 2).
# Here, two morphisms go from one object to another, e.g., f, g: A -> B.
# There are no non-identity morphisms to compose, so composition rules are trivial.
# The structure is unique regardless of how we label the morphisms. This gives 1 category.
count_pattern2 = 1

# Pattern 3: Distribution is (2, 2, 0, 0).
# Hom(A,A) and Hom(B,B) both have 2 morphisms, forming two monoids of order 2.
# There are 2 non-isomorphic monoids of order 2. Let's call them M1 and M2.
# We need to count the number of distinct pairs we can form, up to swapping the objects.
# The possible pairs of monoids are (M1, M1), (M2, M2), and (M1, M2).
# The pair (M2, M1) is isomorphic to (M1, M2). So there are 3 possibilities.
count_pattern3 = 3

# Pattern 4: Distribution is (2, 1, 1, 0) or (1, 2, 0, 1).
# One object has a monoid of order 2 (Hom(A,A)), and there's one morphism from A to B (g: A -> B).
# Let f be the non-identity morphism in Hom(A,A). The composition g . f must equal g.
# The 2 possible monoid structures for Hom(A,A) each yield a valid and distinct category.
count_pattern4 = 2

# Pattern 5: Distribution is (2, 1, 0, 1) or (1, 2, 1, 0).
# Similar to the above, but the morphism goes from B to A (g: B -> A).
# The composition f . g must equal g.
# This structure is not isomorphic to Pattern 4.
# Again, the 2 monoid structures for Hom(A,A) give 2 distinct categories.
count_pattern5 = 2

# Pattern 6: Distribution is (1, 1, 1, 1).
# One morphism f: A -> B, and one morphism g: B -> A.
# The composition rules are uniquely determined by the axioms: g . f = id_A and f . g = id_B.
# This defines a single, rigid category structure where the objects A and B are isomorphic.
count_pattern6 = 1

# The total number of non-isomorphic categories is the sum of counts from all patterns.
total_categories = count_pattern1 + count_pattern2 + count_pattern3 + count_pattern4 + count_pattern5 + count_pattern6

# The final equation shows the contribution from each distinct structural pattern.
print(f"The total number of non-isomorphic categories is the sum of counts from each pattern:")
print(f"{count_pattern1} + {count_pattern2} + {count_pattern3} + {count_pattern4} + {count_pattern5} + {count_pattern6} = {total_categories}")
