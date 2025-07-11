import math

def solve_av_333_3_1324():
    """
    Calculates the number of 1324-avoiding permutations of length 333 with 3 inversions.
    The method is to classify all permutations with 3 inversions, check each class
    for the 1324-avoidance property, count the valid permutations in each class,
    and sum the results.
    """
    n = 333
    k = 3

    print(f"To find av_{n}^{k}(1324), we count the number of 1324-avoiding permutations of length {n} with {k} inversions.")
    print("We classify permutations with 3 inversions into five disjoint types and count the number of avoiding permutations in each type.")
    
    # Type 1: A single block of 3 elements permuted as 'cba'.
    # For example, (..., 3, 2, 1, ...).
    # These contain the 1324 pattern unless they are at the very beginning or end of the permutation.
    # The avoiding permutations correspond to i=1 and i=n-2.
    # e.g., for i=1: (3,2,1,4,5,...,n) is 132-avoiding and thus 1324-avoiding.
    # e.g., for i=n-2: (1,...,n-3,n,n-1,n-2) is a skew-sum of avoider and thus 1324-avoiding.
    count1 = 2
    print(f"\nType 1 ('cba' pattern): Found {count1} avoiding permutations.")

    # Type 2: A single block of 4 elements permuted as 'dabc'.
    # For example, (..., 4, 1, 2, 3, ...).
    # Similar to Type 1, these generally contain 1324 except at the boundaries i=1 and i=n-3.
    count2 = 2
    print(f"Type 2 ('dabc' pattern): Found {count2} avoiding permutations.")

    # Type 3: Two disjoint blocks, 'cab' (2 inversions) and 'ba' (1 inversion).
    # For example, (..., 3, 1, 2, ..., 5, 4, ...).
    # These permutations are decomposable as a direct sum of 1324-avoiding permutations,
    # so they are all 1324-avoiding.
    # We need to choose the starting position 'i' of the 3-element block and 'j' of the 2-element block.
    # Constraints: 1 <= i and i + 2 < j and j <= n-1.
    # The number of ways is Sum_{i=1}^{n-4} (n-1 - (i+3) + 1) = Sum_{i=1}^{n-4} (n-i-3).
    # This sum equals (n-4)*(n-3)/2.
    count3 = (n - 4) * (n - 3) // 2
    print(f"Type 3 ('cab' + 'ba'): Found (n-4)*(n-3)/2 = ({n-4})*({n-3})/2 = {count3} avoiding permutations.")

    # Type 4: Two disjoint blocks, 'bca' (2 inversions) and 'ba' (1 inversion).
    # For example, (..., 2, 3, 1, ..., 5, 4, ...).
    # Like Type 3, these are all 1324-avoiding.
    # The count is the same as for Type 3.
    count4 = (n - 4) * (n - 3) // 2
    print(f"Type 4 ('bca' + 'ba'): Found (n-4)*(n-3)/2 = ({n-4})*({n-3})/2 = {count4} avoiding permutations.")

    # Type 5: Three disjoint 'ba' blocks (1 inversion each).
    # For example, (..., 2, 1, ..., 4, 3, ..., 6, 5, ...).
    # These are also decomposable and thus all 1324-avoiding.
    # We need to choose 3 disjoint pairs of adjacent positions to swap.
    # This is equivalent to choosing 3 indices i,j,k from {1, ..., n-1} such that i+1<j and j+1<k.
    # This corresponds to choosing 3 numbers from {1, ..., n-3}.
    # The count is C(n-3, 3).
    count5 = math.comb(n - 3, 3)
    print(f"Type 5 ('ba' + 'ba' + 'ba'): Found C(n-3, 3) = C({n-3}, 3) = {count5} avoiding permutations.")

    # The five types are disjoint based on the number of non-fixed points in the permutations.
    # Total count is the sum of the counts from each type.
    total_count = count1 + count2 + count3 + count4 + count5

    print("\nSince these five types are disjoint, the total number is their sum:")
    print(f"{count1} + {count2} + {count3} + {count4} + {count5} = {total_count}")

solve_av_333_3_1324()
<<<6043734>>>