import math

def calculate_hamiltonicity_threshold(n, eta):
    """
    Calculates the d-threshold for Hamiltonicity based on n and eta.

    Args:
        n (int): The number of vertices in the graph.
        eta (float): The parameter defining the minimum degree d.

    The problem defines the minimum degree d as:
    d = n/2 - eta
    for 1/2 <= eta <= n/64.

    The d-threshold for Hamiltonicity in this regime is given by the formula:
    p = eta / n
    """

    # --- Parameter validation ---
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (0.5 <= eta <= n / 64):
        print(f"Error: eta = {eta} is not in the valid range [0.5, {n/64:.2f}] for n = {n}.")
        return

    # --- Calculation ---
    # Calculate the minimum degree d
    d = n / 2 - eta

    # The d-threshold probability p
    p_threshold = eta / n

    # --- Output ---
    print(f"Given parameters:")
    print(f"  Number of vertices, n = {n}")
    print(f"  Parameter, eta = {eta}")
    print("-" * 30)
    print("Step 1: Calculate the minimum degree d.")
    # Show the numbers in the equation for d
    # Use math.isclose to handle potential floating point inaccuracies for integer check
    if math.isclose(d, int(d)):
        print(f"  d = n/2 - eta = {n}/2 - {eta} = {int(d)}")
    else:
        print(f"  d = n/2 - eta = {n}/2 - {eta} = {d}")

    print("\nStep 2: Calculate the d-threshold for Hamiltonicity, p.")
    print("  The formula for the threshold is p = eta / n.")
    # Show the numbers in the equation for p
    print(f"  p = {eta} / {n} = {p_threshold}")
    print("-" * 30)
    print("Conclusion:")
    print("The d-threshold for Hamiltonicity is p = eta / n.")
    print(f"For n={n} and eta={eta}, this probability is {p_threshold:.4f}.")


if __name__ == "__main__":
    # Example values satisfying the conditions
    # n=2000, eta=25. Check: 0.5 <= 25 <= 2000/64 = 31.25. Condition holds.
    example_n = 2000
    example_eta = 25
    calculate_hamiltonicity_threshold(example_n, example_eta)
