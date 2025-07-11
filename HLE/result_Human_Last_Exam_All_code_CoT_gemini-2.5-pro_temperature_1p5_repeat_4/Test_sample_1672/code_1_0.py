import numpy as np

def run_entropy_maximization_simulation():
    """
    Simulates an iterative process to maximize state entropy and justifies the answer.
    """
    # Set parameters for the simulation
    NUM_STATES = 10
    NUM_ITERATIONS = 20
    # The 'alpha' parameter models how strongly the new policy seeks out
    # the novel states from the previous iteration.
    ALPHA = 0.5

    def calculate_entropy(p):
        """Calculates the Shannon entropy of a probability distribution."""
        # Add a small epsilon to prevent log(0) for numerical stability
        epsilon = 1e-12
        return -np.sum(p * np.log(p + epsilon))

    # --- Start of the Simulation ---

    # We model the state distribution p_k induced by policy pi^k.
    # Initialize a highly skewed distribution for k=0.
    p = np.arange(NUM_STATES, 0, -1, dtype=float)
    p = p / np.sum(p)  # Normalize to make it a valid probability distribution

    print("This simulation demonstrates how the policy iteration converges to maximize state entropy.")
    print(f"The policy pi^k is trained with rewards r(s) = -log(p_{{k-1}}(s)).")
    print("This encourages visiting less frequent states, driving the distribution towards uniform.")
    print("-" * 70)

    # Track entropy over iterations
    entropy_history = []

    # Iteratively update the distribution
    for k in range(NUM_ITERATIONS):
        entropy = calculate_entropy(p)
        entropy_history.append(entropy)

        # Print status at key iterations
        if k in [0, 1, 5, 10, NUM_ITERATIONS - 1]:
            print(f"Iteration k={k}:")
            # Using list comprehension for cleaner printing of the distribution
            p_str = ', '.join([f'{x:.3f}' for x in p])
            print(f"  State distribution p_{k}: [{p_str}]")
            print(f"  Entropy H(s) = {entropy:.4f}\n")

        # The update rule simulates the outcome of the policy optimization at step k.
        # It creates a new distribution by mixing the old one with a "target"
        # distribution that favors previously infrequent states.
        
        # Target distribution is inversely proportional to the current one.
        epsilon = 1e-12
        p_inverse = 1.0 / (p + epsilon)
        target_dist = p_inverse / np.sum(p_inverse)

        # Update the distribution towards the target
        p = (1 - ALPHA) * p + ALPHA * target_dist

    # --- Analysis of the Final State ---
    print("-" * 70)
    print("Analysis after convergence (as k -> infinity):")

    final_p = p
    final_entropy = entropy_history[-1]
    # The maximum possible entropy for N states is log(N)
    max_entropy = np.log(NUM_STATES)

    print(f"Final state distribution (at k={NUM_ITERATIONS - 1}):")
    final_p_str = ', '.join([f'{x:.4f}' for x in final_p])
    print(f"p_final = [{final_p_str}]")
    uniform_p_str = ', '.join([f'{1/NUM_STATES:.4f}' for _ in range(NUM_STATES)])
    print(f"This is very close to the uniform distribution: [{uniform_p_str}]")

    print("\nFinal Entropy Calculation:")
    # Showing the structure of the final entropy calculation as requested
    term_values = [pi * np.log(pi + 1e-12) for pi in final_p]
    print(f"H(s) = - ( p_0*log(p_0) + p_1*log(p_1) + ... + p_{NUM_STATES-1}*log(p_{NUM_STATES-1}) )")
    print(f"H(s) = - ( {' + '.join([f'{x:.4f}*log({x:.4f})' for x in final_p])} )")
    print(f"H(s) = - ( {''.join([f'({x:.4f})' for x in term_values]).replace(')(',') + (') } )")
    print(f"H(s) = - ( {np.sum(term_values):.4f} )")
    print(f"Final Entropy = {final_entropy:.4f}")

    print(f"\nTheoretical Maximum Entropy = log({NUM_STATES}) = {max_entropy:.4f}")

    print("\nConclusion: As k increases, the state entropy H(s) increases and converges to the")
    print("maximum possible value. This maximum is achieved as k -> infinity.")

run_entropy_maximization_simulation()