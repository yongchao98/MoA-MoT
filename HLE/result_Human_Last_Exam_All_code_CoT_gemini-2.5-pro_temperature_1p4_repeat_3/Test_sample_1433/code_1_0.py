import math

def solve_hypercube_meeting_time(d):
    """
    Calculates the expected time and variance for two random walks to meet on a d-hypercube.
    They start at opposite corners.
    """
    # The distance between walkers is always even, so meeting is only possible for even d.
    if d % 2 != 0:
        return float('inf'), float('inf')

    # Part 1: Calculate the expected time E_d
    # We solve for Delta_k = E_k - E_{k-2} first.
    # E_k is the expected time starting from distance k.
    deltas = {}
    
    # Base case for the recurrence: k = d
    # At distance d, the distance can only decrease or stay the same.
    # The recurrence p_{d,d-2} * Delta_d = 1 gives this.
    deltas[d] = d / (d - 1)
    
    # Recursively calculate Delta_k for k = d-2, d-4, ..., 2
    for k in range(d - 2, 0, -2):
        delta_k_plus_2 = deltas[k + 2]
        
        # The recurrence is: p_{k,k-2} * Delta_k = 1 + p_{k,k+2} * Delta_{k+2}
        # where p_{k,k-2} = k(k-1)/d^2 and p_{k,k+2} = (d-k)(d-k-1)/d^2
        numerator = d**2 + (d - k) * (d - k - 1) * delta_k_plus_2
        denominator = k * (k - 1)
        deltas[k] = numerator / denominator

    # Now, calculate E_k for all k. E_k is the sum of relevant Deltas.
    # E_k = Delta_2 + Delta_4 + ... + Delta_k
    expectations = {}
    current_E = 0
    for k in range(2, d + 1, 2):
        current_E += deltas[k]
        expectations[k] = current_E
    
    E_d = expectations[d]

    # Part 2: Calculate the variance Var(X_d)
    # We first calculate V_k = E[X_k^2], using D_k = V_k - V_{k-2}.
    Ds = {}
    
    # Base case for the recurrence on D_k: k = d
    # p_{d,d-2}*D_d = 2*E_d - 1
    Ds[d] = (d**2 * (2 * E_d - 1)) / (d * (d - 1))

    # Recursively calculate D_k for k = d-2, d-4, ..., 2
    for k in range(d - 2, 0, -2):
        D_k_plus_2 = Ds[k + 2]
        E_k = expectations[k]
        
        # The recurrence is: p_{k,k-2}*D_k = 2*E_k - 1 + p_{k,k+2}*D_{k+2}
        numerator = d**2 * (2 * E_k - 1) + (d - k) * (d - k - 1) * D_k_plus_2
        denominator = k * (k - 1)
        Ds[k] = numerator / den
        
    # V_d = E[X_d^2] is the sum of all D_k
    V_d = sum(Ds.values())
    
    variance_d = V_d - E_d**2
    
    return E_d, variance_d

def check_inequality(d, E_d_val):
    """
    Checks if EX_d <= (d/2) * d^d / d! for a given even d.
    """
    if d % 2 != 0:
        return False
    rhs = (d / 2) * (d**d / math.factorial(d))
    return E_d_val <= rhs

# --- Main Execution ---

# Calculate for d=14
E14, Var14 = solve_hypercube_meeting_time(14)
print(f"EX_14 = {int(E14)}")
print(f"D^2X_14 = {int(Var14)}")

# Answer for d=15
print(f"EX_15 = inf")

# Check the inequality
inequality_holds = True
# We can check for a few even d values to be confident
for d_check in range(2, 16, 2):
    E_d, _ = solve_hypercube_meeting_time(d_check)
    if not check_inequality(d_check, E_d):
        inequality_holds = False
        break
inequality_answer = "yes" if inequality_holds else "no"
print(f"Is it true that for even d number EX_d <= d/2 * d^d/d! ? {inequality_answer}")
