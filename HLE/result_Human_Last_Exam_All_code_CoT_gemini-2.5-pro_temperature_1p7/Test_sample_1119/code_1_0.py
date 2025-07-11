def count_sequences(N, K, M):
    """
    Calculates the number of sequences of K positive integers satisfying the given conditions.

    Args:
      N: The maximum value for any number in the sequence.
      K: The length of the sequence.
      M: The maximum allowed increase between consecutive numbers.

    Returns:
      The total number of possible sequences.
    """
    if M * (K - 1) >= N:
        print(f"Warning: The condition M(K-1) < N (which is {M}*({K}-1) < {N}) is not met.")
        print("The formula might still work, but the problem's premise is violated.")
    
    # We use two arrays to store prefix sums for space optimization.
    # S_prev[j] will store the number of valid sequences of length (i-1) ending with a value <= j.
    # S_curr[j] will store the number of valid sequences of length i ending with a value <= j.
    # We use 1-based indexing for numbers (1 to N), so arrays are size N+1.
    
    # Base case: For sequences of length i=1
    # dp[1][j] = 1 for any j from 1 to N.
    # The prefix sum S_prev[j] = sum_{k=1 to j} dp[1][k] = j.
    S_prev = [j for j in range(N + 1)]
    
    # Iterate from sequence length i = 2 to K
    for i in range(2, K + 1):
        S_curr = [0] * (N + 1)
        for j in range(1, N + 1):
            # Calculate dp[i][j], the number of sequences of length i ending with j.
            # dp[i][j] = sum of dp[i-1][k] for k in [max(1, j-M), j-1]
            # This sum is S_prev[j-1] - S_prev[max(0, j-M-1)]
            
            upper_k = j - 1
            lower_k = max(1, j - M)
            
            dp_i_j = 0
            if lower_k <= upper_k:
                dp_i_j = S_prev[upper_k] - S_prev[lower_k - 1]
            
            # Update the current prefix sum array
            # S_curr[j] = S_curr[j-1] + dp[i][j]
            S_curr[j] = S_curr[j - 1] + dp_i_j
            
        # The current row becomes the previous row for the next iteration
        S_prev = S_curr

    # After the loop, S_prev holds the prefix sums for sequences of length K.
    # The total number of sequences is S_prev[N].
    total_sequences = S_prev[N]
    
    # Reconstruct the last row of dp values (dp[K][j]) from the prefix sums
    dp_K = [0] * (N + 1)
    sum_terms = []
    for j in range(1, N + 1):
        dp_K[j] = S_prev[j] - S_prev[j-1]
        if dp_K[j] > 0:
            sum_terms.append(str(dp_K[j]))

    # Print the detailed breakdown of the final count
    print(f"For N={N}, K={K}, M={M}:")
    print("The number of sequences of length K is the sum of sequences ending with each possible number 'j'.")
    print(f"Contribution of each 'j':")
    for j in range(1, N+1):
      if dp_K[j] > 0:
          print(f"  - Number of sequences ending with {j}: {dp_K[j]}")

    print("\nThe final calculation is:")
    print(f"{' + '.join(sum_terms)} = {total_sequences}")
    
    return total_sequences

if __name__ == '__main__':
    # You can change these values to test with different parameters
    N = 10
    K = 4
    M = 3
    
    # The condition M(K-1) < N is 3*(4-1) = 9 < 10, which holds.
    
    final_count = count_sequences(N, K, M)
    print(f"\nTotal number of possible sequences: {final_count}")


    # Example from the problem description analysis: N=5, K=3, M=2
    # The condition M(K-1) < N is 2*(3-1) = 4 < 5, which holds.
    print("\n" + "="*40)
    print("Another example:")
    count_sequences(N=5, K=3, M=2)
