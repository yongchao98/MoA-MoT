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

def solve_circuit_problem():
    """
    Finds the optimal resistor values based on the problem's constraints
    and calculates the maximum possible current through R3.
    """
    R2 = 6
    
    # We need to find the smallest possible R1 and R3 that satisfy all conditions,
    # as this will maximize the current I3.
    # We will loop through prime numbers for R3 starting from the smallest possible value.
    # From the problem constraints, it can be deduced that R2=6 must be the minimum value.
    
    r3_candidate = R2 + 1
    while True:
        if is_prime(r3_candidate):
            R3 = r3_candidate
            
            # Per the constraint R3 - R1 > 2, R1 must be less than R3 - 2.
            # Also, R2=6 must be the minimum value, so R1 must be greater than 6.
            for R1 in range(R2 + 1, R3 - 2):
                
                # Check the z(C)=R2 condition. This means R2 is the value
                # furthest from the mean of the three resistors.
                mean = (R1 + R2 + R3) / 3.0
                dist_r1 = abs(R1 - mean)
                dist_r2 = abs(R2 - mean)
                dist_r3 = abs(R3 - mean)

                if dist_r2 >= dist_r1 and dist_r2 >= dist_r3:
                    # This is the first valid set of resistors found. Since we are searching
                    # from the smallest possible values, this set gives the maximum current.
                    
                    # Calculate the current I3 using the derived formula.
                    sum_of_products = R1 * R2 + R1 * R3 + R2 * R3
                    numerator = 26 * (R1 + R3) * R2
                    denominator = R3 * sum_of_products
                    
                    # Simplify the fraction to its lowest terms.
                    common_divisor = math.gcd(numerator, denominator)
                    final_num = numerator // common_divisor
                    final_den = denominator // common_divisor

                    print(f"Found optimal resistor values: R1 = {R1} Ω, R2 = {R2} Ω, R3 = {R3} Ω.")
                    print("The maximum possible current through R3 is calculated using the formula:")
                    print("I_3 = (26 * (R1 + R3) * R2) / (R3 * (R1*R2 + R1*R3 + R2*R3))")
                    print("\nSubstituting the values:")
                    print(f"I_3 = (26 * ({R1} + {R3}) * {R2}) / ({R3} * ({R1}*{R2} + {R1}*{R3} + {R2}*{R3}))")
                    print(f"I_3 = (26 * {R1 + R3} * {R2}) / ({R3} * ({R1*R2} + {R1*R3} + {R2*R3}))")
                    print(f"I_3 = ({26 * (R1 + R3) * R2}) / ({R3} * {sum_of_products})")
                    print(f"I_3 = {numerator} / {denominator}")
                    print("\nIn simplest fractional form:")
                    print(f"I_3 = {final_num} / {final_den} Amps")
                    
                    return

        r3_candidate += 1

# Execute the solver function
solve_circuit_problem()