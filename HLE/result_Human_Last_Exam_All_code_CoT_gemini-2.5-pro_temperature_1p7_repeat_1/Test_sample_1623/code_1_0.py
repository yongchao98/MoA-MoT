import itertools

def count_inversions(perm):
    """Counts the number of inversions in a permutation."""
    n = len(perm)
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                inversions += 1
    return inversions

def solve_tb():
    """
    Calculates the Thurston-Bennequin number for the given 5x5 grid diagram.
    """
    n = 5
    
    # 1-based coordinates from the problem statement
    o_coords = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_coords = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]

    # Convert to 0-based permutations tau and sigma
    # tau[i] = row of 'o' in column i
    # sigma[i] = row of 'x' in column i
    tau = [0] * n
    sigma = [0] * n

    for col, row in o_coords:
        tau[col - 1] = row - 1
    
    for col, row in x_coords:
        sigma[col - 1] = row - 1

    print("The grid diagram is defined by two permutations, tau (for 'o') and sigma (for 'x').")
    print(f"tau = {tau}")
    print(f"sigma = {sigma}")
    print("-" * 30)

    # Calculate inversions
    inv_tau = count_inversions(tau)
    inv_sigma = count_inversions(sigma)
    
    print("The Thurston-Bennequin number (tb) is calculated as inv(tau) - inv(sigma).")
    print(f"Number of inversions in tau: {inv_tau}")
    print(f"Number of inversions in sigma: {inv_sigma}")
    print("-" * 30)
    
    # Calculate tb
    tb = inv_tau - inv_sigma

    print("Final Calculation:")
    print(f"tb = {inv_tau} - {inv_sigma} = {tb}")
    
solve_tb()