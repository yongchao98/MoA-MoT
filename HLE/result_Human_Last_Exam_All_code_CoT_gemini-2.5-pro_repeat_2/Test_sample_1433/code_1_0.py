import math
import fractions

def calculate_ex_dx(d):
    """
    Calculates the expected meeting time and its variance on a d-dimensional hypercube.
    d must be an even integer. Returns results as Fraction objects.
    """
    if d % 2 != 0:
        return float('inf'), float('inf')

    # Part 1: Calculate E_k for k = 2, 4, ..., d
    # The state is the Hamming distance k. Let E_k be the expected time to reach state 0 from k.
    # The recurrence for the differences Delta_k = E_k - E_{k-2} is derived as:
    # Delta_{k+2} = (k^2 * Delta_k - d^2) / (d-k)^2
    num_states = d // 2
    # delta[j] corresponds to Delta_{2j}
    delta = [fractions.Fraction(0)] * (num_states + 1)
    # Ek[j] corresponds to E_{2j}
    Ek = [fractions.Fraction(0)] * (num_states + 1)
    
    # Boundary condition at k=d: E_d = 1 + E_{d-2}, so Delta_d = 1
    delta[num_states] = fractions.Fraction(1)
    
    # Solve recurrence backwards for delta_j (where j=k/2)
    # delta_j = ((d-2j)^2 * delta_{j+1} + d^2) / (2j)^2
    d_sq = d**2
    for j in range(num_states - 1, 0, -1):
        k = 2 * j
        numerator = (d - k)**2 * delta[j+1] + d_sq
        denominator = k**2
        delta[j] = numerator / denominator
        
    # E_d is the sum of all Delta_k
    Ek[0] = fractions.Fraction(0) # E_0 = 0
    for j in range(1, num_states + 1):
        Ek[j] = Ek[j-1] + delta[j]
        
    exd = Ek[num_states]

    # Part 2: Calculate S_k = E[X^2] for k = 2, 4, ..., d
    # The recurrence for Psi_k = S_k - S_{k-2} is:
    # Psi_{k+2} = (k^2 * Psi_k - d^2 * (2*E_k - 1)) / (d-k)^2
    psi = [fractions.Fraction(0)] * (num_states + 1)
    
    # Boundary condition at k=d: S_d = S_{d-2} + 2*E_d - 1, so Psi_d = 2*E_d - 1
    psi[num_states] = 2 * Ek[num_states] - 1
    
    # Solve recurrence backwards for psi_j (where j=k/2)
    # psi_j = ((d-2j)^2 * psi_{j+1} + d^2*(2*E_{2j}-1)) / (2j)^2
    for j in range(num_states - 1, 0, -1):
        k = 2 * j
        numerator = (d - k)**2 * psi[j+1] + d_sq * (2 * Ek[j] - 1)
        denominator = k**2
        psi[j] = numerator / denominator
        
    # S_d is the sum of all Psi_k
    Sk = [fractions.Fraction(0)] * (num_states + 1) # S_0 = 0
    for j in range(1, num_states + 1):
        Sk[j] = Sk[j-1] + psi[j]
        
    exd_sq = Sk[num_states]
    varxd = exd_sq - exd**2
    
    return exd, varxd

# --- Main execution ---

# Task 1 & 2: Calculate EX_14 and D^2X_14
d14 = 14
Ex14_frac, VarX14_frac = calculate_ex_dx(d14)
Ex14_int = int(float(Ex14_frac))
VarX14_int = int(float(VarX14_frac))

print(f"The integer part of the expected time to meet for d=14 is:")
print(f"EX_14 = {Ex14_int}")

print(f"\nThe integer part of the variance for d=14 is:")
print(f"D^2X_14 = {VarX14_int}")

# Task 3: Calculate EX_15
print(f"\nThe expected time to meet for d=15 is:")
print("EX_15 = inf")

# Task 4: Verify the inequality for d=14
# EX_d <= (d/2) * (d^d / d!)
d = 14
lhs = float(Ex14_frac)
# Use math.loggamma for precision with large factorials
# log(d!) = lgamma(d+1)
log_rhs = math.log(d / 2.0) + d * math.log(d) - math.lgamma(d + 1)
rhs = math.exp(log_rhs)

is_true = lhs <= rhs
answer_yes_no = "yes" if is_true else "no"

print(f"\nIs it true that for even d, EX_d <= (d/2) * (d^d / d!)?")
print(f"We check for d={d}:")
print(f"The equation is EX_14 <= (14/2) * (14^14 / 14!)")
print(f"The integer part of the left side (EX_14) is {int(lhs)}")
print(f"The integer part of the right side is {int(rhs)}")
print(f"Is {int(lhs)} <= {int(rhs)}? {answer_yes_no.upper()}")
