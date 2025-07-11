import itertools

def solve():
    """
    Calculates the maximal Thurston-Bennequin number for a given grid diagram.
    """
    n = 5
    # O positions: (1,1), (2,2), (3,3), (4,4), (5,5)
    # X positions: (1,4), (2,5), (3,1), (4,2), (5,3)

    # 1. Represent X's as a permutation X_i = column of X in row i
    # X in row 1 is at (3,1), so X_1 = 3
    # X in row 2 is at (4,2), so X_2 = 4
    # X in row 3 is at (5,3), so X_3 = 5
    # X in row 4 is at (1,4), so X_4 = 1
    # X in row 5 is at (2,5), so X_5 = 2
    X = [3, 4, 5, 1, 2]
    print(f"The grid size is n = {n}")
    print(f"The permutation for X markers is: {X}")
    print("-" * 20)

    # 2. Calculate the writhe (w)
    # inv(X) is the number of inversions in the permutation X
    inversions = 0
    inv_pairs = []
    for i, j in itertools.combinations(range(n), 2):
        if X[i] > X[j]:
            inversions += 1
            inv_pairs.append((X[i], X[j]))
            
    print(f"Step 1: Calculate the number of inversions in X.")
    print(f"The inversion pairs are: {inv_pairs}")
    print(f"Number of inversions, inv(X) = {inversions}")
    
    total_pairs = n * (n - 1) // 2
    writhe = 2 * inversions - total_pairs
    print(f"\nStep 2: Calculate the writhe w(G).")
    print(f"The formula for writhe is w(G) = 2 * inv(X) - n(n-1)/2")
    print(f"w(G) = 2 * {inversions} - {total_pairs} = {writhe}")
    print("-" * 20)

    # 3. Calculate the Thurston-Bennequin number for the grid (tb(G))
    tb_G = writhe - n
    print(f"Step 3: Calculate the Thurston-Bennequin number tb(G) for this grid.")
    print(f"The formula is tb(G) = w(G) - n")
    print(f"tb(G) = {writhe} - {n} = {tb_G}")
    print("-" * 20)

    # 4. Count the number of destabilizations
    num_destabilizations = 0
    destab_indices = []
    for i in range(n - 1):
        if X[i+1] == X[i] + 1:
            num_destabilizations += 1
            destab_indices.append(i+1) # Using 1-based indexing for output
            
    print(f"Step 4: Count the number of possible destabilizations.")
    print(f"A destabilization exists for each index i (1 to {n-1}) where X[i+1] = X[i] + 1.")
    print(f"Destabilizations found at indices i = {destab_indices}")
    print(f"Total number of destabilizations = {num_destabilizations}")
    print("-" * 20)

    # 5. Calculate the maximal Thurston-Bennequin number
    tb_max = tb_G + num_destabilizations
    print(f"Step 5: Calculate the maximal Thurston-Bennequin number TB(K).")
    print(f"The formula is TB(K) = tb(G) + (number of destabilizations)")
    print(f"The final equation is:")
    print(f"{tb_max} = {tb_G} + {num_destabilizations}")
    
    # Return the final answer
    return tb_max

if __name__ == '__main__':
    result = solve()
    print(f"\nThe maximal Thurston-Bennequin number is {result}.")
    print(f"<<<{result}>>>")
