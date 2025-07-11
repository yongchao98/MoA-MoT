import numpy as np

def simulate_polynomial_cot():
    """
    Simulates a transformer with polynomial steps of Chain-of-Thought (CoT).

    This function illustrates why such a process is in the complexity class P.
    The simulation runs for a number of steps polynomial to the input size 'n'.
    """
    # n: Represents the initial input size (e.g., sequence length)
    n = 10

    # d: Represents the model's embedding dimension
    d = 32

    # Let the number of CoT steps be polynomial in 'n'. We'll use n for simplicity.
    polynomial_steps = n

    print(f"Starting CoT simulation with:")
    print(f"  Input size n = {n}")
    print(f"  Embedding dimension d = {d}")
    print(f"  Polynomial steps = {polynomial_steps}\n")

    # A fixed matrix representing the 'reasoning' computation of one transformer step
    # The size is d x (n+d) to handle a growing context
    # In a real scenario, this would be more complex (attention, etc.)
    # For this simulation, we'll use a simpler d x d matrix for clarity
    transformer_matrix = np.random.rand(d, d)

    # Initial 'thought' or context vector (e.g., from the input prompt)
    # We'll use a vector of size 'd' for simplicity
    current_thought = np.random.rand(d)
    
    total_operations = 0

    # The core loop: runs a polynomial number of times
    for i in range(polynomial_steps):
        # This simulates one step of reasoning.
        # The computation here is a matrix-vector product, which takes O(d*d) time.
        # In a real transformer, it would be poly(n, d).
        previous_thought = current_thought
        current_thought = np.dot(transformer_matrix, previous_thought)

        # The key aspect of CoT: the output is used in the next step.
        # Here, 'current_thought' becomes the input for the next iteration.
        
        operations_in_step = transformer_matrix.shape[0] * transformer_matrix.shape[1]
        total_operations += operations_in_step

        print(f"Step {i + 1}/{polynomial_steps}: Generated new 'thought' vector.")
        # print(current_thought) # Uncomment to see the vector

    print("\n--- Simulation Complete ---")
    print(f"The simulation performed {polynomial_steps} sequential steps of computation.")
    print("Since the number of steps is polynomial in 'n' (it is 'n' here),")
    print("and each step takes polynomial time, the entire process is in P.")

    # The user requested to print the numbers in a final equation.
    # We define the total operations as Steps * (Operations per Step).
    op_per_step_str = f"(d*d = {d*d})"
    num_steps = polynomial_steps
    print("\nFinal Complexity Equation:")
    print(f"Total Operations = (Number of Steps) * (Operations per Step)")
    print(f"Total Operations = {num_steps} * {op_per_step_str} = {total_operations}")


if __name__ == '__main__':
    simulate_polynomial_cot()