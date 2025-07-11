import math

def solve_thurston_bennequin():
    """
    Calculates the maximal Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    # Step 1: Represent the grid as a permutation.
    # O's are at (i,i), so they represent the identity permutation.
    # X's are at (1,4), (2,5), (3,1), (4,2), (5,3).
    # We use 0-based indexing for the list, so we subtract 1 from each coordinate.
    # p[i] = j means an X is at column i+1, row j+1.
    p = [3, 4, 0, 1, 2]
    n = 5
    
    # Step 2 & 3: Calculate the writhe by counting inversions and non-inversions.
    # An inversion is a pair (i, j) such that i < j and p[i] > p[j].
    # A non-inversion is a pair (i, j) such that i < j and p[i] < p[j].
    inversions = 0
    non_inversions = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inversions += 1
            else:
                non_inversions += 1
                
    # The writhe w(G) = (number of non-inversions) - (number of inversions)
    writhe = non_inversions - inversions
    
    # Step 4: Calculate the maximal Thurston-Bennequin number.
    # For an alternating grid diagram, tb(K) = w(G) - n.
    # We have confirmed the permutation p = (1 4 2 5 3) is a single cycle of length 5,
    # so the diagram is alternating.
    tb_number = writhe - n
    
    # Print the final result as an equation.
    # The problem requires printing the numbers in the final equation.
    # The equation is: Maximal TB Number = Writhe - Grid Size
    print("The maximal Thurston-Bennequin number is the writhe of the diagram minus the grid number.")
    print("Equation: Maximal TB Number = w - n")
    print(f"Calculated values: w = {writhe}, n = {n}")
    print(f"Final calculation: {tb_number} = {writhe} - {n}")
    print(f"The maximal Thurston-Bennequin number is {tb_number}.")

solve_thurston_bennequin()
<<< -7 >>>