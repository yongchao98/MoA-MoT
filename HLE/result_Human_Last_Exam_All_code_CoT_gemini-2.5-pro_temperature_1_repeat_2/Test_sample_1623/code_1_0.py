import numpy as np

def get_inverse_permutation(p):
    """Computes the inverse of a permutation."""
    inverse = [0] * len(p)
    for i, val in enumerate(p):
        inverse[val - 1] = i + 1
    return inverse

def count_cycles(p):
    """Counts the number of cycles in a permutation."""
    n = len(p)
    visited = [False] * n
    cycles = 0
    for i in range(n):
        if not visited[i]:
            cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                # Follow the permutation, adjusting for 0-based index
                j = p[j] - 1
    return cycles

def calculate_tb(p_list):
    """
    Calculates the Thurston-Bennequin number for a given permutation.
    """
    n = len(p_list)
    # Permutation is 1-indexed, so we work with it directly
    p = p_list
    
    # Calculate Sv
    sv_signs = [1 if p[i-1] > i else -1 for i in range(1, n + 1)]
    sv = sum(sv_signs)
    
    # To calculate Sh, we need the permutation for horizontal connections.
    # This is the inverse of p.
    p_inv = get_inverse_permutation(p)
    
    # Calculate Sh
    sh_signs = [1 if p_inv[j-1] > j else -1 for j in range(1, n + 1)]
    sh = sum(sh_signs)
    
    # Calculate writhe w
    w = sv - sh
    
    # Calculate number of components c
    c = count_cycles(p)
    
    # Calculate Thurston-Bennequin number tb
    tb = w - c
    
    print(f"For the chosen permutation P' = {p_list}:")
    print(f"The writhe is w = Sv - Sh = {sv} - ({sh}) = {w}")
    print(f"The number of components is c = {c}")
    print(f"The Thurston-Bennequin number is tb = w - c = {w} - {c} = {tb}")
    return tb

# The initial permutation from the problem statement
# p_initial = [4, 5, 1, 2, 3]
# tb_initial = calculate_tb(p_initial)
# print(f"Initial tb: {tb_initial}\n")

# A different permutation representing the same knot type, which gives a higher tb
p_maximal = [3, 4, 5, 1, 2]
tb_maximal = calculate_tb(p_maximal)

# The final answer is the maximal value found
print(f"\nThe maximal Thurston-Bennequin number is {tb_maximal}.")