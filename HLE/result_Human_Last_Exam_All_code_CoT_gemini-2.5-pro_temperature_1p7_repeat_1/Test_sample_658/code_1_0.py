import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value of d_{B,delta}.
    
    Args:
        B (float): The upper bound of the interval [0, B].
        delta (float): The maximum allowed error.
        
    Returns:
        float: The asymptotic value of the minimum degree d.
    """
    if B <= 0 or delta <= 0 or delta >= 1:
        print("Error: B must be > 0, and delta must be in (0, 1).")
        return None

    # L is a common notation in approximation theory for log(1/delta)
    L = math.log(1/delta)
    
    # The asymptotic formula is composed of two terms that dominate in different regimes.
    # Term 1: Dominates when B >> L
    term1 = math.sqrt(B * L)
    
    # Term 2: Dominates when L >> B. The argument of the log is B/L + L/B
    # to keep it well-defined and capture the correct asymptotic behavior.
    log_argument = (B / L) + (L / B)
    term2 = L / math.log(log_argument)
    
    asymptotic_d = term1 + term2
    
    return asymptotic_d, L, log_argument

# --- User-defined values ---
B_val = 100.0
delta_val = 0.001
# ---------------------------

result = calculate_asymptotic_degree(B_val, delta_val)

if result is not None:
    d_approx, L_val, log_arg_val = result
    
    print("This script calculates the asymptotic value for d_{B,delta}.")
    print("The formula used is: d ≈ sqrt(B*L) + L / log(B/L + L/B), where L = log(1/delta).\n")
    print(f"For the input values B = {B_val} and delta = {delta_val}:")
    print(f"L = log(1/{delta_val}) = {L_val}")
    
    print("\nThe equation with these values plugged in is:")
    print(f"d ≈ sqrt({B_val} * {L_val}) + {L_val} / log({B_val}/{L_val} + {L_val}/{B_val})")
    print(f"d ≈ sqrt({B_val * L_val}) + {L_val} / log({log_arg_val})")
    
    print("\nEvaluated result:")
    print(f"d ≈ {d_approx}")
