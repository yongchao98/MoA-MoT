def solve():
    """
    Calculates av_{333}^3(1324), the number of 1324-avoiding permutations
    of length 333 with 3 inversions.
    """
    n = 333

    # Case 1: Permutations with a ...cba... structure on 3 consecutive elements.
    # These are 132-avoiding. Number of choices for the block is n-2.
    count1 = n - 2

    # Case 2: Permutations with a ...dabc... structure on 4 consecutive elements.
    # These are 132-avoiding. Number of choices for the block is n-3.
    count2 = n - 3
    
    # Case 3: Permutations with a ...bcda... structure on 4 consecutive elements.
    # These are 132-avoiding. Number of choices for the block is n-3.
    count3 = n - 3

    # Case 4: Permutations with a ...bdac... structure on 4 consecutive elements.
    # These contain 132 and are 1324-avoiding only if the elements are {n-3,...,n}.
    # This gives 1 permutation.
    count4 = 1

    # Case 5: Permutations with a ...cadb... structure on 4 consecutive elements.
    # These contain 132 and are 1324-avoiding only if the elements are {n-3,...,n}.
    # This gives 1 permutation.
    count5 = 1

    # Case 6: Permutations with two disjoint blocks of inversions (2+1).
    # These contain 132 and are 1324-avoiding only if the 1-inversion block is {n-1,n}.
    # The 2-inversion block on {i,i+1,i+2} can be chosen in n-4 ways.
    count6 = n - 4

    total_count = count1 + count2 + count3 + count4 + count5 + count6

    print(f"{count1} + {count2} + {count3} + {count4} + {count5} + {count6} = {total_count}")

solve()