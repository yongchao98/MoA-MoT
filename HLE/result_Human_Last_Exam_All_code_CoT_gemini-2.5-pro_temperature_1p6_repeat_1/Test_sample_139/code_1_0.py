import math

def is_prime(n):
    """Checks if a number n is prime."""
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

# Initialize variables to store the solution
R1_sol, R3_sol = None, None
# Given constants
R2 = 6
V_fail = 26

# To maximize the current, we need the smallest possible resistor values.
# We derived the conditions: R1 > 8 and R1 + 2 < R3 < 2*R1 - 6.
# We will search for the smallest R1 that yields a valid prime R3.
for R1 in range(9, 100):  # Search for R1 starting from 9
    r3_lower_bound = R1 + 2
    r3_upper_bound = 2 * R1 - 6

    # Search for a prime R3 in the calculated range
    for R3 in range(r3_lower_bound + 1, r3_upper_bound):
        if is_prime(R3):
            R1_sol = R1
            R3_sol = R3
            break  # Found the smallest valid R3 for this R1
    if R1_sol is not None:
        break  # Found the smallest R1 with a solution, so stop

if R1_sol and R3_sol:
    # Calculate the total source current from the failure condition
    # I_total = V_fail / R_eq_fail = V_fail / ((R1 * R3) / (R1 + R3))
    I_total = V_fail * (R1_sol + R3_sol) / (R1_sol * R3_sol)

    # Calculate the current through R3 using the current divider rule for the intact circuit
    # I_R3 = I_total * G3 / (G1 + G2 + G3), where G = 1/R
    G1 = 1 / R1_sol
    G2 = 1 / R2
    G3 = 1 / R3_sol
    current_r3 = I_total * G3 / (G1 + G2 + G3)

    print("The optimal resistor values that satisfy the conditions are:")
    print(f"R1 = {R1_sol} Ohms")
    print(f"R2 = {R2} Ohms")
    print(f"R3 = {R3_sol} Ohms")
    print("\nThe equation for the current through R3 (I_R3) with these values is:")
    print(f"I_R3 = ({V_fail} * ({R1_sol} + {R3_sol}) / ({R1_sol} * {R3_sol})) * (1/{R3_sol}) / (1/{R1_sol} + 1/{R2} + 1/{R3_sol})")
    
    print("\nThe maximum possible current through R3 is:")
    print(f"{current_r3}")

else:
    print("No valid resistor values found within the search range.")