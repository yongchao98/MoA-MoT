import math

def calculate_asymptotic_d(B, delta):
    """
    Calculates and prints the asymptotic value A(B, delta) for d_{B,delta}.
    
    The asymptotic formula is A(B, delta) = B + L^2 / B, where L = log(1/delta).
    This function demonstrates the formula for given B and delta values.
    """
    if B < 1 or delta <= 0 or delta >= 1:
        print(f"Invalid input for B={B}, delta={delta}. B must be >= 1 and delta must be in (0, 1).")
        return

    L = math.log(1 / delta)
    term_B = B
    term_L2_B = (L**2) / B
    A = term_B + term_L2_B

    print(f"For B = {B} and delta = {delta}:")
    print(f"We calculate L = log(1/delta) = log({1/delta:.2e}) = {L:.4f}")
    print("The asymptotic value is given by the formula: B + L^2 / B")
    # The following line prints each number in the final equation
    print(f"Result: {term_B:.4f} + ({L:.4f}^2 / {B:.4f}) = {term_B:.4f} + {term_L2_B:.4f} = {A:.4f}")
    print("-" * 30)

# --- Examples ---

# Example 1: B is large compared to L (A is dominated by the B term)
print("Regime 1: B is large, L is small")
calculate_asymptotic_d(1000, 0.01)

# Example 2: L is large compared to B (A is dominated by the L^2/B term)
print("Regime 2: L is large, B is small")
calculate_asymptotic_d(2, 1e-20)

# Example 3: Crossover regime where both terms are of similar magnitude (B is comparable to L)
print("Regime 3: Crossover (B and L are comparable)")
# We choose B=10 and L=10, which corresponds to delta = e^-10
calculate_asymptotic_d(10, math.exp(-10))
