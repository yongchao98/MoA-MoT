import math
from fractions import Fraction

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

# This part of the code performs a search for the resistor values
# based on the logical constraints derived from the problem statement.
R2 = 6
found_solution = False
r1_sol, r3_sol = 0, 0

# We need to find the smallest valid (R1, R3) pair to maximize the current.
# The condition (R3 + 6)/2 < R1 < R3 - 2 implies R3 must be a prime > 10.
# We will search for the smallest such prime R3, starting from 11.
r3_candidate = 11
while not found_solution:
    if is_prime(r3_candidate):
        # Calculate bounds for R1
        r1_lower_bound = (r3_candidate + 6) / 2
        r1_upper_bound = r3_candidate - 2

        # Check if an integer R1 exists in the interval (lower, upper).
        # We find the smallest integer R1 in that range.
        r1_candidate = math.floor(r1_lower_bound) + 1
        if r1_candidate < r1_upper_bound:
            # Found the smallest (R1, R3) pair satisfying all constraints.
            r1_sol = r1_candidate
            r3_sol = r3_candidate
            found_solution = True
            
    if not found_solution:
        # Move to the next potential prime candidate.
        r3_candidate += 1

# --- Calculation and Output Section ---

# Calculate the source current I_s from the failure condition.
# When R2 fails open, V_fail = 26 V across the parallel R1 and R3.
V_fail = 26
Is_frac = Fraction(V_fail * (r1_sol + r3_sol), (r1_sol * r3_sol))

# Calculate the current I3 through R3 when R2 is intact using the current divider rule.
G1 = Fraction(1, r1_sol)
G2 = Fraction(1, R2)
G3 = Fraction(1, r3_sol)
Total_Conductance = G1 + G2 + G3
I3_frac = Is_frac * G3 / Total_Conductance

print("Step-by-step Derivation of the Solution:")
print("\n1. Finding Resistor Values:")
print(f"   - Based on the problem's constraints, we derived the inequality (R3 + 6)/2 < R1 < R3 - 2, where R3 is a prime > 10.")
print(f"   - To maximize the final current, we must find the smallest possible resistance values.")
print(f"   - The smallest prime R3 > 10 is 13. For R3 = 13, the inequality becomes 9.5 < R1 < 11, which gives a unique integer solution R1 = 10.")
print(f"   - The determined resistor values are: R1 = {r1_sol} Ω, R2 = {R2} Ω, R3 = {r3_sol} Ω.")

print("\n2. Calculating Maximum Possible Current through R3:")
print("   - Assuming a parallel circuit, we first find the source current I_s from the failure condition:")
print(f"     I_s * (R1 * R3) / (R1 + R3) = 26")
print(f"     I_s * ({r1_sol} * {r3_sol}) / ({r1_sol} + {r3_sol}) = {V_fail}")
print(f"     I_s = {V_fail} * ({r1_sol} + {r3_sol}) / ({r1_sol} * {r3_sol}) = {Is_frac.numerator}/{Is_frac.denominator} A ≈ {float(Is_frac):.3f} A")
print("\n   - With I_s known, we calculate the current I3 through R3 (when R2 is intact) using the current divider rule:")
print("     I3 = I_s * (1/R3) / (1/R1 + 1/R2 + 1/R3)")
print(f"     I3 = ({Is_frac.numerator}/{Is_frac.denominator}) * (1/{r3_sol}) / (1/{r1_sol} + 1/{R2} + 1/{r3_sol})")
print(f"   - Final Calculation:")
print(f"     I3 = ({Is_frac.numerator}/{Is_frac.denominator}) A * ({G3.numerator}/{G3.denominator} S) / ({Total_Conductance.numerator}/{Total_Conductance.denominator} S) = {I3_frac.numerator}/{I3_frac.denominator} A")

print(f"\nThe maximum possible current through R3 is {I3_frac.numerator}/{I3_frac.denominator} A, which is approximately {float(I3_frac):.4f} A.")
