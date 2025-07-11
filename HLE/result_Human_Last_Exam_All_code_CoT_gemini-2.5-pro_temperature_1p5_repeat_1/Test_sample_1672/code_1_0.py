import numpy as np

def simulate_entropy_maximization():
    """
    Simulates the iterative process of policy updates based on an intrinsic
    reward that encourages visiting less frequent states, and shows that this
    process maximizes state entropy.
    """
    # 1. Setup the environment
    N_STATES = 5
    # Initial policy leads to a skewed state distribution
    # State 0 is visited much more often than others.
    p = np.array([0.8, 0.05, 0.05, 0.05, 0.05])
    
    ITERATIONS = 20
    # Alpha is a learning rate or mixing factor to ensure smooth convergence
    ALPHA = 0.5
    # Epsilon to prevent log(0)
    EPSILON = 1e-9

    # 2. Calculate the maximum possible entropy for a uniform distribution
    uniform_p = np.full(N_STATES, 1.0 / N_STATES)
    max_entropy = -np.sum(uniform_p * np.log(uniform_p))

    print(f"Number of states: {N_STATES}")
    print(f"Maximum possible entropy (for uniform distribution): {max_entropy:.4f}\n")
    print("Starting simulation...")
    
    # 3. Run the iterative simulation
    for k in range(ITERATIONS):
        # Calculate the entropy of the current state distribution
        entropy = -np.sum(p * np.log(p + EPSILON))
        
        # Output the state of the system at the current iteration
        # The instruction "output each number in the final equation" is interpreted
        # as showing the numerical values at each step of the simulation.
        print(f"Iteration {k:2d}: p = {np.round(p, 4)}, Entropy = {entropy:.4f}")

        # The reward for visiting a state is -log(p(s))
        # The target distribution q(s) is proportional to exp(reward), which is 1/p(s)
        # This is the distribution that maximally explores states based on the last policy.
        q_unnormalized = 1.0 / (p + EPSILON)
        q = q_unnormalized / np.sum(q_unnormalized)

        # The new policy's state distribution is a mix of the old and the target.
        # This represents a learning step towards the exploratory target.
        p = (1 - ALPHA) * p + ALPHA * q

    # Final state after all iterations
    final_entropy = -np.sum(p * np.log(p + EPSILON))
    print(f"Iteration {ITERATIONS:2d}: p = {np.round(p, 4)}, Entropy = {final_entropy:.4f}\n")

    print("Conclusion from simulation:")
    print("As the number of iterations (k) increases, the state distribution 'p' becomes")
    print("more uniform, and its entropy approaches the maximum possible value.")
    print("This demonstrates that the limiting policy (as k -> infinity) is the one that")
    print("maximizes the state entropy.")


if __name__ == '__main__':
    simulate_entropy_maximization()
