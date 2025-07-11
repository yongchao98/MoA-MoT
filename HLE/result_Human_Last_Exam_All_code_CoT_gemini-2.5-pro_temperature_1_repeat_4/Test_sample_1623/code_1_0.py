import numpy as np

def solve():
    """
    Calculates the Thurston-Bennequin number for a knot from a given grid diagram.
    """
    n = 5
    # The positions are given in (column, row) format, 1-indexed.
    O_pos = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    X_pos = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # Convert to 0-indexed permutations for programming convenience.
    # sigma_O[i] = j means an 'o' in column i+1 is in row j+1.
    sigma_O = [p[1] - 1 for p in sorted(O_pos)]
    sigma_X = [p[1] - 1 for p in sorted(X_pos)]

    # The inverse permutation sigma_O_inv is needed.
    # Since O's are on the main diagonal, sigma_O is the identity, as is its inverse.
    sigma_O_inv = list(range(n))

    # The inverse permutation sigma_X_inv is also needed.
    sigma_X_inv = [0] * n
    for i in range(n):
        sigma_X_inv[sigma_X[i]] = i

    # Create a set of all marked positions for quick lookup.
    marked_pos = set((p[0] - 1, p[1] - 1) for p in O_pos + X_pos)

    # 1. Calculate the writhe (w)
    w = 0
    for c in range(n):  # column index
        for r in range(n):  # row index
            if (c, r) not in marked_pos:
                v_sign = np.sign(sigma_X[c] - sigma_O[c])
                h_sign = np.sign(sigma_X_inv[r] - sigma_O_inv[r])
                term = v_sign * h_sign
                w += term

    # 2. Calculate the number of right cusps (#R)
    # This is the number of local maxima in the sigma_O permutation.
    num_R = 0
    if n == 1:
        num_R = 1
    else:
        # Check first point
        if sigma_O[0] > sigma_O[1]:
            num_R += 1
        # Check middle points
        for i in range(1, n - 1):
            if sigma_O[i - 1] < sigma_O[i] and sigma_O[i] > sigma_O[i + 1]:
                num_R += 1
        # Check last point
        if sigma_O[n - 1] > sigma_O[n - 2]:
            num_R += 1
            
    # 3. Calculate the Thurston-Bennequin number (tb)
    tb = w - num_R

    # Print the results step-by-step
    print(f"The calculation for the Thurston-Bennequin number is tb = w - #R.")
    print(f"1. The writhe w is calculated by summing the signs of all {n*n - 2*n} crossings.")
    print(f"   The calculated writhe is w = {w}.")
    print(f"2. The number of right cusps #R is the number of local maxima in the 'o' permutation σ_O = {[x+1 for x in sigma_O]}.")
    print(f"   For σ_O, the only local maximum is at the last position, so #R = {num_R}.")
    print(f"3. The final Thurston-Bennequin number is:")
    print(f"   tb = w - #R = {w} - {num_R} = {tb}")

solve()