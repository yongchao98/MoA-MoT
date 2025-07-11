import numpy as np

def calculate_chain_quality(beta, p, N=50):
    """
    Calculates the expected chain quality in a Bitcoin system with a selfish miner.

    Args:
        beta (float): The fraction of mining power controlled by the adversary (0 < beta < 1).
        p (float): The probability of choosing the adversary's block in a tie-breaking
                   situation when an honest miner finds the next block (0 <= p <= 1).
        N (int): The maximum lead state to consider for the numerical solution.

    Returns:
        float: The expected chain quality.
    """
    if not (0 < beta < 1):
        raise ValueError("beta must be between 0 and 1")
    if not (0 <= p <= 1):
        raise ValueError("p must be between 0 and 1")
        
    alpha = 1.0 - beta

    # The system of linear equations is Ax = b, where x is the vector of unknowns.
    # x = [h_1, h_2, ..., h_N, l_1, l_2, ..., l_N]
    # So the total number of variables is 2N.
    A = np.zeros((2 * N, 2 * N))
    b = np.zeros(2 * N)

    # --- Equations for h_i (expected honest blocks) ---
    # h_1 = beta * h_2 + alpha * h_0_fork
    h_0_fork = alpha * (p * 1.0 + (1.0 - p) * 2.0)
    A[0, 0] = 1.0
    A[0, 1] = -beta
    b[0] = alpha * h_0_fork

    # h_2 = beta * h_3 + alpha * 0  (Adversary publishes 2 blocks, 0 honest)
    A[1, 1] = 1.0
    if N > 2:
      A[1, 2] = -beta
    b[1] = 0

    # h_k = beta * h_{k+1} + alpha * h_{k-2} for k >= 3
    # This comes from the assumption that the attacker's lead reduces by 2.
    for k in range(3, N + 1):
        idx = k - 1
        A[idx, idx] = 1.0
        if k < N:
            A[idx, idx + 1] = -beta
        # h_{k-2} index is k-3. h_0 is not in our system, but h_1 is index 0
        A[idx, idx - 2] = -alpha
        b[idx] = 0.0

    # --- Equations for l_i (expected chain length increase) ---
    # l_1 = beta * l_2 + alpha * l_0_fork
    # l_0_fork is always 2.
    l_0_fork = 2.0
    offset = N
    A[offset, offset] = 1.0
    A[offset, offset + 1] = -beta
    b[offset] = alpha * l_0_fork
    
    # l_2 = beta * l_3 + alpha * 2 (Adversary publishes 2 blocks)
    A[offset + 1, offset + 1] = 1.0
    if N > 2:
        A[offset + 1, offset + 2] = -beta
    b[offset + 1] = alpha * 2.0

    # l_k = beta * l_{k+1} + alpha * (1 + l_{k-2}) for k >= 3
    # This comes from the assumption that the attacker publishes 1 block and lead -> k-2.
    for k in range(3, N + 1):
        idx = offset + k - 1
        A[idx, idx] = 1.0
        if k < N:
            A[idx, idx + 1] = -beta
        A[idx, idx - 2] = -alpha
        b[idx] = alpha * 1.0

    # Solve the system
    try:
        x = np.linalg.solve(A, b)
        h_1 = x[0]
        l_1 = x[N]
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations. The matrix might be singular.")
        return None

    # Calculate final expected chain quality
    # E[H] = alpha * 1 + beta * h_1
    # E[L] = alpha * 1 + beta * l_1
    expected_honest_blocks = alpha + beta * h_1
    expected_total_length = alpha + beta * l_1

    if expected_total_length == 0:
        return 0.0

    chain_quality = expected_honest_blocks / expected_total_length
    
    print(f"For beta = {beta} and p = {p}:")
    print(f"The expected number of honest blocks per round, E[H] = {expected_honest_blocks:.4f}")
    print(f"The expected total blocks per round, E[L] = {expected_total_length:.4f}")
    print("\nThe expected chain quality is the ratio of these two values.")
    print("\nFinal Calculation:")
    print(f"Expected Chain Quality = ({alpha} + {beta} * {h_1:.4f}) / ({alpha} + {beta} * {l_1:.4f})")
    print(f"                       = {expected_honest_blocks:.4f} / {expected_total_length:.4f}")
    print(f"                       = {chain_quality:.4f}")

    return chain_quality

if __name__ == '__main__':
    # Example Parameters
    beta = 0.35  # Adversary's mining power
    p = 0.5     # Tie-breaking probability

    # --- Run the calculation ---
    quality = calculate_chain_quality(beta, p)
    if quality is not None:
        print(f"\n<<<The final expected chain quality is {quality:.4f}>>>")
