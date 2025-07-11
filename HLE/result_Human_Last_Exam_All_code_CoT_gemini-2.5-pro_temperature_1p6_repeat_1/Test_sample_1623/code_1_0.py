import collections

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a knot from a grid diagram.
    """
    # 1. Define grid parameters based on the problem description.
    n = 5
    
    # 2. Define the permutation pi from the X coordinates.
    # The X's are at (1,4), (2,5), (3,1), (4,2), (5,3).
    # This corresponds to the 1-based permutation [4, 5, 1, 2, 3].
    # We convert it to a 0-based list for Python indexing.
    # pi[column_index] = row_index
    pi = [3, 4, 0, 1, 2]
    
    # 3. Calculate l(G).
    # l(G) is the number of pairs (i, j) such that i < j and pi[i] < pi[j].
    l_G = 0
    for i in range(n):
        for j in range(i + 1, n):
            if pi[i] < pi[j]:
                l_G += 1
                
    # 4. Calculate the maximal Thurston-Bennequin number using the formula tb_max = l(G) - n.
    tb_max = l_G - n
    
    # 5. Output the explanation and the final calculation.
    print("The formula for the maximal Thurston-Bennequin number (tb_max) is:")
    print("tb_max = l(G) - n")
    print(f"\nGrid size n = {n}")
    pi_1_based = [x + 1 for x in pi]
    print(f"Permutation π from X's (1-based) = {pi_1_based}")
    print(f"Calculated value l(G) = {l_G}")
    
    print("\nFinal calculation:")
    print(f"{l_G} - {n} = {tb_max}")
    print(f"\nThe maximal Thurston-Bennequin number is {tb_max}.")

solve()