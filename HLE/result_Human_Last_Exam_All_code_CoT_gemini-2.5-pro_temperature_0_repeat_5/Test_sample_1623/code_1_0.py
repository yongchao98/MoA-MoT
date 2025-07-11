import numpy as np

def calculate_num_components(perm, n):
    """
    Calculates the number of cycles in a permutation.
    """
    visited = [False] * n
    cycles = 0
    for i in range(n):
        if not visited[i]:
            cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = perm[j]
    return cycles

def calculate_writhe(sigma, pi, n):
    """
    Calculates the writhe of the grid diagram.
    """
    # Calculate inverse permutations
    sigma_inv = [0] * n
    pi_inv = [0] * n
    for i in range(n):
        sigma_inv[sigma[i]] = i
        pi_inv[pi[i]] = i

    writhe = 0
    print("Calculating writhe by summing signs of all crossings:")
    # Iterate over all grid cells (i, j)
    for i in range(n):
        for j in range(n):
            # A crossing exists if the cell (i, j) is empty
            if j != sigma[i] and j != pi[i]:
                # sgn(vertical orientation) * sgn(horizontal orientation)
                sign_v = 1 if pi[i] > sigma[i] else -1
                sign_h = 1 if pi_inv[j] > sigma_inv[j] else -1
                crossing_sign = sign_v * sign_h
                writhe += crossing_sign
    return writhe

def solve_grid_diagram_tb():
    """
    Solves for the Thurston-Bennequin number of the given 5x5 grid diagram.
    """
    n = 5
    # o's at (1,1), (2,2), (3,3), (4,4), (5,5)
    # Using 0-based indexing:
    sigma = [0, 1, 2, 3, 4]
    
    # x's at (1,4), (2,5), (3,1), (4,2), (5,3)
    # Using 0-based indexing:
    pi = [3, 4, 0, 1, 2]

    print(f"Grid size n = {n}")
    print(f"Permutation for o's (sigma): {sigma}")
    print(f"Permutation for x's (pi):    {pi}")
    print("-" * 30)

    # Step 1: Calculate number of components (c)
    # c = number of cycles in pi * sigma_inv
    # Since sigma is identity, sigma_inv is identity, so we just need cycles in pi.
    print("Step 1: Calculate the number of components (c)")
    c = calculate_num_components(pi, n)
    print(f"The number of components c is the number of cycles in permutation pi, which is {c}.")
    print("-" * 30)

    # Step 2: Calculate writhe (w)
    print("Step 2: Calculate the writhe (w)")
    w = calculate_writhe(sigma, pi, n)
    print(f"The total writhe w is {w}.")
    print("-" * 30)

    # Step 3: Calculate Thurston-Bennequin number (tb)
    print("Step 3: Calculate the Thurston-Bennequin number (tb = w - c)")
    tb = w - c
    print(f"The maximal Thurston-Bennequin number is the writhe minus the number of components.")
    print(f"tb = w - c = {w} - {c} = {tb}")
    
    return tb

if __name__ == '__main__':
    final_answer = solve_grid_diagram_tb()
    print(f"\nFinal Answer: {final_answer}")
    # The final answer is wrapped in <<<>>> as requested.
    # print(f"<<<{final_answer}>>>")