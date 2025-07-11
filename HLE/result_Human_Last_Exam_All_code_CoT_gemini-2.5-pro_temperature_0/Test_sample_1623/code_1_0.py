import sys

def solve_knot_tb():
    """
    Calculates the maximal Thurston-Bennequin number for the knot
    associated with the given grid diagram.
    """
    # The grid size is 5x5.
    n = 5
    
    # Step 1: Define the permutations from the grid positions.
    # The problem uses 1-based indexing (columns and rows from 1 to 5).
    # We will convert to 0-based indexing for calculations.
    # O's are at (1,1), (2,2), (3,3), (4,4), (5,5).
    # This means the permutation pi_O, where pi_O[col-1] = row-1, is the identity.
    # pi_O = [0, 1, 2, 3, 4]
    
    # X's are at (1,4), (2,5), (3,1), (4,2), (5,3).
    # This gives the permutation pi_X.
    # pi_X[0] = 3, pi_X[1] = 4, pi_X[2] = 0, pi_X[3] = 1, pi_X[4] = 2
    pi_X = [3, 4, 0, 1, 2]
    
    # The knot is determined by the permutation pi = pi_X * pi_O_inverse.
    # Since pi_O is the identity, its inverse is also the identity.
    # So, the knot permutation is simply pi_X.
    pi = pi_X
    
    # Step 2: Identify the knot type.
    # A grid diagram with pi_O as identity and pi_X(i) = (i + k) mod n
    # represents the torus knot T(n, k).
    # Let's find k for our permutation.
    # For i=0, pi[0] = 3. So, (0 + k) mod 5 = 3, which means k=3.
    # We can verify this for other indices:
    # (1+3)%5 = 4 (correct)
    # (2+3)%5 = 0 (correct)
    # (3+3)%5 = 1 (correct)
    # (4+3)%5 = 2 (correct)
    # So, the knot is the T(5, 3) torus knot.
    p = 5
    q = 3
    
    print(f"The grid diagram corresponds to the T({p},{q}) torus knot.")
    
    # Step 3: Use the formula for the maximal Thurston-Bennequin number.
    # The formula for a T(p,q) torus knot (with p,q > 0) is tb_max = p*q - p - q.
    print("The maximal Thurston-Bennequin number for a T(p,q) torus knot is given by the formula: tb_max = p*q - p - q.")
    
    # Step 4: Calculate the result.
    tb_max = p * q - p - q
    
    print("\nFor the T({},{}) knot, the calculation is:".format(p, q))
    print("tb_max = {} * {} - {} - {}".format(p, q, p, q))
    print("tb_max = {} - {}".format(p * q, p + q))
    print("tb_max = {}".format(tb_max))
    
    # The final answer is tb_max.
    # To conform to the output format, we will print the final value separately.
    # This is handled outside the function call.
    return tb_max

if __name__ == '__main__':
    solve_knot_tb()