import scipy.integrate
import sys

# Check if scipy is installed, as it's required for the calculation.
try:
    import scipy.integrate
except ImportError:
    print("Error: The 'scipy' library is required to run this script.")
    print("Please install it using: pip install scipy")
    sys.exit(1)

def solve_probability():
    """
    Explains and calculates the probability of an edge in a jointly exchangeable random graph.
    """
    # --- Step 1: Explain the Theoretical Framework ---
    print("--- Theoretical Framework ---")
    print("For a random graph that is jointly exchangeable, its structure is described by the Aldous-Hoover representation theorem.")
    print("This theorem states that such a graph can be generated via a latent variable model:")
    print("1. Each node 'i' is assigned a hidden attribute alpha_i, drawn independently from a Uniform distribution on [0, 1]. This corresponds to the 'uniform measure U'.")
    print("2. The probability of an edge y_ij between nodes i and j is determined by a function W(alpha_i, alpha_j). This function W, mapping [0,1]^2 to [0,1], is called a graphon.")
    print("3. The term 'random measure F' implies that the graphon W can itself be random, drawn from a distribution over possible graph-generating functions.")
    print("\n")

    # --- Step 2: Explain the Probability Calculation ---
    print("--- Calculating the Probability ---")
    print("The unconditional probability of an edge, P(y_ij = 1), is the same for all pairs (i, j) due to exchangeability.")
    print("It is found by averaging the edge probability W(u, v) over all random choices.")
    print("This is done by integrating the expected graphon, E[W(u, v)], over the unit square where the latent variables live.")
    print("The general formula is:")
    print("P(y_ij = 1) = Integral from v=0 to 1 of [Integral from u=0 to 1 of E[W(u, v)] du] dv")
    print("\n")

    # --- Step 3: Demonstrate with a Concrete Example ---
    print("--- Example Calculation ---")
    print("Since the specific graphon W is not provided, we will demonstrate the calculation with an example.")
    print("Let's assume the expected graphon has a simple form: E[W(u, v)] = (u + v) / 2.0")
    print("\n")

    # Define the example expected graphon function
    def expected_graphon(u, v):
      """
      An example of an expected graphon function, E[W(u, v)].
      """
      return (u + v) / 2.0

    # Define integration limits, which are the 'numbers in the equation'
    lower_bound = 0.0
    upper_bound = 1.0

    # Perform the double integration using scipy.integrate.dblquad
    # The function returns a tuple (result, error_estimate)
    try:
        probability, _ = scipy.integrate.dblquad(
            expected_graphon,
            lower_bound,
            upper_bound,
            lambda u: lower_bound, # Lower limit for the inner integral (v)
            lambda u: upper_bound  # Upper limit for the inner integral (v)
        )
    except Exception as e:
        print(f"An error occurred during integration: {e}")
        return

    # --- Step 4: Output the Final Equation and Result ---
    print("--- Final Equation and Result ---")
    print("For our example, the equation we are solving is:")
    # We print each number that defines the equation
    print(f"P(y_ij=1) = Integral from v={lower_bound} to {upper_bound} [Integral from u={lower_bound} to {upper_bound} [(u + v) / 2.0 du] dv]")
    print("\n")
    print("The calculated probability for this example is:")
    print(f"P(y_ij = 1) = {probability}")

if __name__ == "__main__":
    solve_probability()