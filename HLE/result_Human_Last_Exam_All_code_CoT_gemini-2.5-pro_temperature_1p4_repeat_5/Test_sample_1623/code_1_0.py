import numpy as np

def calculate_thurston_bennequin_number():
    """
    Calculates the maximal Thurston-Bennequin number for a knot specified by a grid diagram.

    The problem provides a 5x5 grid with 'o's and 'x's at specific positions.
    The Thurston-Bennequin number (tb) for the Legendrian knot representative
    from this grid diagram is given by the formula:
        tb = trace(pi) - n + 2 * inversions(pi)
    where:
    - n is the grid size.
    - pi is the knot permutation derived from the grid.
    - trace(pi) is the number of fixed points of the permutation.
    - inversions(pi) is the number of inversions in the permutation.
    """
    # Grid size
    n = 5

    # Positions of 'o's (1-indexed from problem)
    # (1,1), (2,2), (3,3), (4,4), (5,5)
    # Converted to 0-indexed (col, row)
    O_pos = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)]

    # Positions of 'x's (1-indexed from problem)
    # (1,4), (2,5), (3,1), (4,2), (5,3)
    # Converted to 0-indexed (col, row)
    X_pos = [(0, 3), (1, 4), (2, 0), (3, 1), (4, 2)]

    # Step 1: Construct the permutations sigma_O and sigma_X
    # sigma_O maps column i to the row of 'o' in that column
    # sigma_X maps column i to the row of 'x' in that column
    sigma_O = {col: row for col, row in O_pos}
    sigma_X = {col: row for col, row in X_pos}

    # Step 2: Calculate the inverse permutation sigma_O_inv
    # sigma_O_inv maps a row j to the column of 'o' in that row
    sigma_O_inv = {row: col for col, row in sigma_O.items()}
    
    # Step 3: Compute the knot permutation pi = sigma_X * sigma_O_inv
    # pi[i] = sigma_X(sigma_O_inv(i))
    pi_list = [sigma_X[sigma_O_inv[i]] for i in range(n)]
    
    print(f"The knot permutation pi is: {pi_list}")

    # Step 4: Calculate the trace of pi (number of fixed points)
    trace = 0
    for i in range(n):
        if pi_list[i] == i:
            trace += 1
    
    print(f"The trace of the permutation (number of fixed points) is: {trace}")

    # Step 5: Calculate the number of inversions of pi
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if pi_list[i] > pi_list[j]:
                inversions += 1
    
    print(f"The number of inversions in the permutation is: {inversions}")
    print(f"The grid size n is: {n}")

    # Step 6: Apply the formula to get the Thurston-Bennequin number
    tb = trace - n + 2 * inversions
    
    print("\nThe Thurston-Bennequin number is calculated using the formula:")
    print("tb = trace(pi) - n + 2 * inversions(pi)")
    print("\nPlugging in the values, we get:")
    print(f"tb = {trace} - {n} + 2 * {inversions} = {tb}")

calculate_thurston_bennequin_number()
<<<7>>>