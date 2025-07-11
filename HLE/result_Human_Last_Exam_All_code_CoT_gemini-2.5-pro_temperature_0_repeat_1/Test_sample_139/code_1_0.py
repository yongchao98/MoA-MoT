import math

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

def find_resistors_and_current():
    """
    Finds the values of R1 and R3 that satisfy the problem conditions
    and calculates the maximum possible current through R3.
    """
    R2 = 6
    
    # We need to find the smallest R1 and R3 that satisfy the conditions to maximize I3.
    # Our analysis shows R1, R3 > 6 and R1 > 8. We start the search from R1 = 9.
    for R1 in range(9, 100):  # Search over a reasonable range for R1
        # R3 must be a prime number and satisfy R3 - R1 > 2, so R3 >= R1 + 3.
        # The z(C) analysis also gives R3 < 2*R1 - 6.
        for R3 in range(R1 + 3, 2 * R1 - 6):
            if not is_prime(R3):
                continue

            # Check the z(C)=6 condition, which means 6 is the value furthest from the mean.
            # This translates to the following inequalities:
            # |6 - mean| > |R1 - mean|  and  |6 - mean| > |R3 - mean|
            # which simplifies to:
            # |12 - R1 - R3| > |2*R1 - R3 - 6|  and  |12 - R1 - R3| > |2*R3 - R1 - 6|
            
            dist_mean_6 = abs(12 - R1 - R3)
            dist_mean_R1 = abs(2 * R1 - R3 - 6)
            dist_mean_R3 = abs(2 * R3 - R1 - 6)

            if dist_mean_6 > dist_mean_R1 and dist_mean_6 > dist_mean_R3:
                # This is the first valid pair (R1, R3) found. Since we are searching
                # from the smallest R1 upwards, this pair will yield the maximum current.
                
                # Calculate the current I3 using the derived formula:
                # I3 = 156 * (R1 + R3) / (R1*R3 + 6*R1 + 6*R3)
                numerator = 156 * (R1 + R3)
                denominator = (R1 * R3 + 6 * R1 + 6 * R3)
                max_current_I3 = numerator / denominator

                print("Found the resistor values that satisfy all conditions:")
                print(f"R1 = {R1} ohms")
                print(f"R2 = {R2} ohms")
                print(f"R3 = {R3} ohms (which is a prime number)")
                print("\nThese values result in the maximum possible current through R3.")
                print("\nThe equation for the current I3 is:")
                print(f"I3 = 156 * (R1 + R3) / (R1*R3 + 6*R1 + 6*R3)")
                print("\nSubstituting the found values into the equation:")
                print(f"I3 = 156 * ({R1} + {R3}) / ({R1}*{R3} + 6*{R1} + 6*{R3})")
                print(f"I3 = {numerator} / {denominator}")
                print(f"\nThe maximum possible current through R3 is {max_current_I3:.3f} A.")
                
                return max_current_I3

    return None

# Run the solver and capture the final answer
final_answer = find_resistors_and_current()
if final_answer:
    print(f"<<<{final_answer:.3f}>>>")
