import numpy as np
from scipy.special import psi

def calculate_expected_watermark_score(prob_distributions):
    """
    Calculates the exact expected watermarking score E[S] for a sequence 
    of token probability distributions from a language model.

    The formula for the score at a single position t is:
    E_t = sum_{i=1 to K} [ p_{t,i} * (gamma + psi(1 + 1/p_{t,i})) ]
    where gamma is the Euler-Mascheroni constant and psi is the digamma function.
    The total score E[S] is the sum of E_t over all positions t.

    Args:
        prob_distributions (list of list of floats): A list where each element 
                                                     is a probability distribution D_t 
                                                     for a token. Example: [[0.9, 0.1], [0.5, 0.5]]
    
    Returns:
        None. It prints the result.
    """
    # The Euler-Mascheroni constant can be defined as -psi(1)
    euler_gamma = -psi(1)
    
    total_expected_score = 0.0
    total_entropy = 0.0
    n = len(prob_distributions)

    if n == 0:
        print("The list of probability distributions is empty.")
        return

    print("Calculating token by token...")
    # Header for the table
    print("-" * 65)
    print(f"{'Token (t)':<12}{'Entropy H_t':<15}{'Expected Score E_t':<20}{'Gain (E_t - 1)':<15}")
    print("-" * 65)

    for i, p_dist in enumerate(prob_distributions):
        p_t = np.array(p_dist)
        
        # Ensure probabilities sum to 1 and handle potential floating point inaccuracies
        if not np.isclose(np.sum(p_t), 1.0):
            print(f"Warning: Probabilities for token {i+1} do not sum to 1. Normalizing.")
            p_t = p_t / np.sum(p_t)

        # Filter out probabilities that are zero to avoid division by zero and log(0)
        p_t_nonzero = p_t[p_t > 0]
        
        # Calculate Entropy H_t for the current token
        entropy_t = -np.sum(p_t_nonzero * np.log(p_t_nonzero))
        total_entropy += entropy_t
        
        # Calculate expected score E_t for the current token
        e_t = np.sum(p_t_nonzero * (euler_gamma + psi(1.0 + 1.0 / p_t_nonzero)))
        total_expected_score += e_t

        # Print the results for the current token
        print(f"{i+1:<12}{entropy_t:<15.4f}{e_t:<20.4f}{e_t - 1:<15.4f}")

    # Calculate average entropy alpha
    alpha = total_entropy / n

    print("-" * 65)
    print("\n--- Summary ---")
    print(f"Number of tokens (n): {n}")
    print(f"Average entropy per token (alpha): {alpha:.4f} nats")
    print(f"Expected score for non-watermarked text (E[S_unwatermarked] = n): {float(n)}")
    print(f"Total Expected Score for Watermarked Text (E[S]): {total_expected_score:.4f}")
    print(f"Total Expected Gain (E[S] - n): {total_expected_score - n:.4f}")


# --- Example Usage ---
# You can replace this with your own data.
# This example has a text with n=4 tokens and a vocabulary of size K=5.
example_prob_distributions = [
    [0.2, 0.2, 0.2, 0.2, 0.2],    # High entropy (uniform)
    [0.6, 0.1, 0.1, 0.1, 0.1],    # Medium entropy
    [0.94, 0.02, 0.02, 0.01, 0.01], # Low entropy
    [0.5, 0.5, 0.0, 0.0, 0.0]     # Binary choice
]

calculate_expected_watermark_score(example_prob_distributions)