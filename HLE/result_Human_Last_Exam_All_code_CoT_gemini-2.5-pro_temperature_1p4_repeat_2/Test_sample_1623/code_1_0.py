import math

def calculate_tb_from_grid():
    """
    Calculates the Thurston-Bennequin number for a knot from a specific grid diagram.
    """
    # Step 1 & 2: Define grid size and permutation from the problem
    # n=5 grid
    # o's at (1,1), (2,2), (3,3), (4,4), (5,5) -> identity permutation
    # x's at (1,4), (2,5), (3,1), (4,2), (5,3)
    # 1-based permutation sigma = (4, 5, 1, 2, 3)
    # Using 0-based indexing for the list in Python
    n = 5
    sigma = [3, 4, 0, 1, 2]

    # Step 3: Calculate the writhe (w)
    # First, find the number of inversions in sigma
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if sigma[i] > sigma[j]:
                inversions += 1
    
    # Then, calculate the writhe
    writhe = (n * (n - 1) // 2) - (2 * inversions)

    # Step 4: Calculate the number of right cusps (#R)
    # First, find the inverse permutation sigma_inv
    sigma_inv = [0] * n
    for i in range(n):
        sigma_inv[sigma[i]] = i
        
    # Count right cusps from o's: sigma[i] < i
    right_cusps_o = sum(1 for i in range(n) if sigma[i] < i)
    
    # Count right cusps from x's: sigma_inv[j] > j
    right_cusps_x = sum(1 for j in range(n) if sigma_inv[j] > j)
    
    num_right_cusps = right_cusps_o + right_cusps_x

    # Step 5: Calculate the Thurston-Bennequin number (tb)
    # The formula is tb = w - #R / 2
    tb = writhe - num_right_cusps / 2
    
    # Output the results in the required format
    print(f"The calculation for the Thurston-Bennequin number (tb) follows the formula: tb = w - #R/2")
    print(f"1. The writhe (w) is calculated as n(n-1)/2 - 2*inv(σ) = {n*(n-1)//2} - 2*{inversions} = {writhe}")
    print(f"2. The number of right cusps (#R) is |{{i|σ(i)<i}}| + |{{j|σ⁻¹(j)>j}}| = {right_cusps_o} + {right_cusps_x} = {num_right_cusps}")
    print(f"3. Plugging these values into the formula:")
    print(f"tb = {writhe} - {num_right_cusps}/2 = {int(tb)}")


calculate_tb_from_grid()
<<< -5 >>>