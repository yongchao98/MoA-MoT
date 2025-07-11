import numpy as np

def solve():
    """
    Calculates the maximum number of 9x1x1 blocks in a 15x15x15 cube
    using a coloring argument.
    """
    edge_length = 15
    block_length = 9

    # The number of unit cubes with a given color can be found using
    # polynomial multiplication.
    # Let P(t) = 1 + t + t^2 + ... + t^14. The coefficient of t^s in P(t)^3
    # is the number of integer solutions to x+y+z = s where 0 <= x,y,z < 15.

    # 1. Define the polynomial for one dimension.
    # The coefficients are all 1 for powers from 0 to 14.
    p_coeffs = np.ones(edge_length, dtype=np.int64)

    # 2. Cube the polynomial to represent three dimensions.
    # This is done by convolving the polynomial with itself twice.
    p2_coeffs = np.convolve(p_coeffs, p_coeffs)
    p3_coeffs = np.convolve(p2_coeffs, p_coeffs)

    # 3. Calculate the color counts N(k).
    # N(k) is the sum of all coefficients of p3_coeffs whose power
    # is congruent to k modulo 9.
    color_counts = np.zeros(block_length, dtype=np.int64)
    for i, count in enumerate(p3_coeffs):
        color_counts[i % block_length] += count

    # 4. The maximum number of blocks is the minimum of the color counts.
    max_blocks = np.min(color_counts)

    print("The number of unit cubes for each color k=(x+y+z)%9 are:")
    for i, count in enumerate(color_counts):
        print(f"Color {i}: {count} cubes")

    print("\nThe maximum number of blocks is the minimum of these counts.")
    
    # Format the final equation string
    counts_str = ", ".join(map(str, color_counts))
    print(f"max_blocks = min({counts_str})")
    print(f"max_blocks = {max_blocks}")

solve()