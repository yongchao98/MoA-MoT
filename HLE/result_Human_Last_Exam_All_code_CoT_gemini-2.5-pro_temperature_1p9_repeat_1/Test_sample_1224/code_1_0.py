import math

def gcd(a, b):
    """Helper function to calculate the greatest common divisor."""
    return math.gcd(a, b)

def titan_mul(frac1, frac2):
    """
    Simulates a multiplication operation on the Titan 4-bit architecture.
    A fraction is represented as a tuple (numerator, denominator).
    It checks if the resulting numerator or denominator exceeds 15.
    """
    n1, d1 = frac1
    n2, d2 = frac2

    # --- Pre-simplification based on rules ---
    # We can simplify across fractions, e.g., (a/b) * (c/d) = (a/d) * (c/b)
    # This is equivalent to dividing out common factors before multiplication.
    g1 = gcd(n1, d2)
    n1, d2 = n1 // g1, d2 // g1

    g2 = gcd(n2, d1)
    n2, d1 = n2 // g2, d1 // g2
    
    # --- Perform Multiplication ---
    res_n = n1 * n2
    res_d = d1 * d2

    # --- Constraint Check ---
    # All numerators and denominators must remain within the 4-bit range (0-15).
    # This check applies to the final result of the single multiplication operation.
    if res_n > 15 or res_d > 15:
        print(f"Failed to multiply {frac1} by {frac2}.")
        print(f"  Intermediate simplified product: ({n1}*{n2})/({d1}*{d2}) = {res_n}/{res_d}")
        raise ValueError(f"Overflow Error! Result {res_n}/{res_d} exceeds 4-bit limit.")
    
    print(f"Successfully multiplied {frac1} by {frac2} to get {res_n}/{res_d}")
    return (res_n, res_d)

def main():
    """
    Main function to attempt the calculation for Pandora's gravity.
    """
    print("--- Titan Computer Feasibility Test ---")
    print("Attempting to calculate a component of Pandora's volume: 4 * r^3\n")

    try:
        # Define constants as 4-bit fractions (mantissas)
        # Using the simplest possible approximations to give Titan the best chance.
        # Sphere volume constant from (4/3)*pi with pi=3 is 4.
        four = (4, 1)
        
        # Pandora's equatorial radius is ~2000 km. We only use the mantissa.
        # r = 2e6 m. Mantissa is 2.
        r_mantissa = (2, 1)

        # We need to calculate 4 * r * r * r
        print(f"Calculating 4 * r^3 where r has a mantissa of {r_mantissa}...")
        
        # First multiplication: 4 * r
        step1_res = titan_mul(four, r_mantissa)

        # Second multiplication: (4*r) * r
        step2_res = titan_mul(step1_res, r_mantissa)

        # Third multiplication: (4*r*r) * r
        step3_res = titan_mul(step2_res, r_mantissa)
        
        print("\nAll calculations succeeded without overflow.")

    except ValueError as e:
        print(f"\nCOMPUTATION FAILED: {e}")
        print("\nConclusion: The Titan architecture cannot perform the calculation.")
        print("Even basic multiplication of the physical constants results in numbers")
        print("larger than what a 4-bit system can represent. Therefore, calculating")
        print("the planet's mass, gravity, and the final fall time is not feasible.")

if __name__ == "__main__":
    main()
