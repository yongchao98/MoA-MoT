import numpy as np

def solve_and_explain():
    """
    Solves for the maximal rank of the Tongan flag matrix and explains the reasoning.
    """

    # 1. Theoretical Analysis
    # Let 'a' be the value for red pixels and 'b' for white pixels.
    # The flag's matrix has two distinct types of rows:
    #
    # - Row Type 1: An all-red row. This occurs for any row fully in the red field
    #   or intersecting the horizontal bar of the red cross in the canton.
    #   It looks like: [a, a, a, ..., a]
    #
    # - Row Type 2: A row passing through the white sections of the canton.
    #   It contains 'b's for the white background and 'a's for the red cross/field.
    #   It looks like: [b, ..., b, a, ..., a, b, ..., b, a, ..., a]
    #
    # The entire row space of the matrix is spanned by these two vectors.
    # Therefore, the rank can be at most 2.

    # To find the *maximal* rank, we must choose 'a' and 'b' such that these
    # two row types are linearly independent.
    #
    # They are linearly dependent if one is a scalar multiple of the other.
    # For this to happen (with non-zero vectors), they would need to be identical
    # in pattern, which would require a = b. If a = b, the flag is one color,
    # and the rank is 1.
    #
    # By choosing a != b (and a != 0), we make the two row types linearly independent.
    # For example, let a = 2 and b = 1.
    # A vector of all 2s cannot be a multiple of a vector containing 1s and 2s.
    # Since the two row types are linearly independent, the rank is 2.

    # 2. Demonstration with NumPy
    # We can create a sample matrix and compute its rank.
    # Let's choose values to maximize the rank.
    a = 2.0
    b = 1.0

    # Define matrix and feature dimensions (proportions can be approximate)
    H, W = 20, 40
    canton_h, canton_w = H // 2, W // 4
    matrix = np.full((H, W), fill_value=a)  # Red field

    # White canton
    matrix[0:canton_h, 0:canton_w] = b

    # Red cross in canton
    cross_h_bar_start, cross_h_bar_end = canton_h // 3, (canton_h * 2) // 3
    cross_v_bar_start, cross_v_bar_end = canton_w // 3, (canton_w * 2) // 3
    matrix[cross_h_bar_start:cross_h_bar_end, 0:canton_w] = a
    matrix[0:canton_h, cross_v_bar_start:cross_v_bar_end] = a

    # Calculate the rank
    maximal_rank = np.linalg.matrix_rank(matrix)

    # 3. Final Answer
    # The final equation is simply that the maximal rank is this calculated value.
    print("The final equation is: Maximal Rank = {}".format(int(maximal_rank)))
    print("The number in the final equation is:")
    print(int(maximal_rank))

solve_and_explain()