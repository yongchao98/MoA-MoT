import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def solve_circuit():
    """
    Solves for the resistor values and calculates the current through R3.
    """
    R2 = 6
    found_solution = None

    # R3 must be a prime factor of 26
    possible_R3 = [p for p in range(1, 27) if 26 % p == 0 and is_prime(p)]

    # Iterate through possible values of R3 and a reasonable range for R1
    for R3 in possible_R3:
        I = 26 / R3
        # R1 is a distinct positive integer. Let's check a reasonable range.
        for R1 in range(1, 51):
            # Condition 1: Resistors are distinct integers
            if R1 == R2 or R1 == R3:
                continue

            # Condition 2: R3 - R1 > 2
            if R3 - R1 <= 2:
                continue

            # Condition 3: z(R1, R2, R3) = R2
            # This means R2 is the value furthest from the mean.
            resistors = [R1, R2, R3]
            mean = sum(resistors) / 3
            dist_R1 = abs(R1 - mean)
            dist_R2 = abs(R2 - mean)
            dist_R3 = abs(R3 - mean)

            # Check if R2's distance to the mean is the greatest
            # Add a small tolerance for floating point comparisons if necessary,
            # but with integers it should be exact.
            if dist_R2 > dist_R1 and dist_R2 > dist_R3:
                found_solution = {'R1': R1, 'R2': R2, 'R3': R3, 'I': I}
                # Since we found a unique solution from our derivation, we can stop.
                break
        if found_solution:
            break
            
    if found_solution:
        # Calculate the current through R3 when R2 is intact using the current divider rule
        I = found_solution['I']
        R2 = found_solution['R2']
        R3 = found_solution['R3']
        
        I3 = I * (R2 / (R2 + R3))
        
        print(f"Found a unique valid solution for the components:")
        print(f"R1 = {found_solution['R1']} ohms")
        print(f"R2 = {R2} ohms")
        print(f"R3 = {R3} ohms")
        print(f"I_source = {I} A\n")

        print(f"The current through R3 when R2 is intact is calculated using the current divider formula:")
        print(f"I3 = I_source * (R2 / (R2 + R3))")
        print(f"I3 = {I} A * ({R2} / ({R2} + {R3}))")
        print(f"I3 = {I} A * ({R2} / {R2 + R3})")
        print(f"I3 = {I * R2} / {R2 + R3}")
        print(f"The maximum possible current through R3 is: {I3:.4f} A")

    else:
        print("No solution found that satisfies all conditions.")

solve_circuit()
<<<0.6316>>>