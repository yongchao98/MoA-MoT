import numpy as np

def calculate_asymptotic_d(B, delta):
    """
    Calculates the asymptotic value of d_{B,delta}.

    Args:
      B: The upper bound of the interval [0, B].
      delta: The desired accuracy.

    Returns:
      The asymptotic value of d.
    """
    if B < 1:
        # The problem states B >= 1
        B = 1
    if not (0 < delta < 1):
        raise ValueError("delta must be in (0, 1)")

    L = np.log(1/delta)
    
    # We add 2 in the log to keep the argument > 1, this doesn't affect the
    # Theta notation for the asymptotic behavior.
    # The asymptotic behavior is a combination of two regimes.
    
    # Behavior for B >> L
    term1 = np.sqrt(B * L)
    
    # Behavior for L >> B. The argument of log must be > 1.
    # We use max(2, L/B) which is a simple way to ensure this.
    log_arg = L / B
    if log_arg < 2:
      log_arg = 2
    term2 = L / np.log(log_arg)
    
    # A robust approximation for the asymptotic value is the sum of the
    # dominant terms from each regime.
    asymptotic_d = term1 + term2
    
    print("The asymptotic value is a function of B and L = log(delta^-1).")
    print("A(B, L) = sqrt(B*L) + L / log(max(2, L/B))")
    print(f"For B = {B} and delta = {delta} (L = {L:.4f}):")
    print(f"sqrt(B*L) = {term1:.4f}")
    print(f"L / log(max(2, L/B)) = {term2:.4f}")
    # The final expression is presented piece by piece for clarity, but it represents a single formula.
    # We combine them into the final answer.
    
    # Let's re-state the combined expression as the result.
    print(f"The final asymptotic value A(B, L) is the sum of these terms, which for these values is approximately: {asymptotic_d:.4f}")

# You can modify these values to test different scenarios
B_val = 100.0
delta_val = 1e-10

calculate_asymptotic_d(B_val, delta_val)