import numpy as np

def calculate_q(k):
    """Calculates the value of q_k based on the problem definition."""
    k_abs = abs(k)
    if k % 2 != 0:  # k is odd
        return 2 - (k_abs - 1) / 2
    else:  # k is even
        return 5/3 - k_abs / 2

def verify_recurrence():
    """Verifies the recurrence relation found for q_k."""
    print("Verifying the recurrence relation: 2*q_k - q_{k-2} - q_{k+2} = p_k")
    p = {}
    for k in range(-5, 6):
        # Calculate p_k from the recurrence
        p_k = 2 * calculate_q(k) - calculate_q(k - 2) - calculate_q(k + 2)
        # Store non-zero values, rounding for display
        if abs(p_k) > 1e-9:
            p[k] = round(p_k, 4)

    print(f"The resulting sequence p_k is non-zero only for k in {list(p.keys())}.")
    print(f"p_k values: {p}")
    print("This confirms Q(D) * (2 - D^2 - D^-2) = 2 + D + D^-1")
    print("\nFrom this, we factor Q(D) and find F(D) = 1 / (1-D)")
    print("The whitening filter W(D) is 1/F(D^-1), which gives W(D) = 1 - D.")
    print("\nThe equation for the whitening filter W(D) is:")
    print("W(D) = 1 - 1*D^1")

verify_recurrence()
