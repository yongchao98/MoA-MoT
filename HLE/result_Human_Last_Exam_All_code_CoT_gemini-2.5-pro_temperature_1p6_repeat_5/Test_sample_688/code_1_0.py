import math

def calculate_prefactor_cn():
    """
    Calculates the system-independent prefactor c_n for the fully f-connected
    Ree-Hoover diagram in the virial expansion.

    The formula for the prefactor is c_n = -(n-1) / n!.
    Since n is not specified in the problem, we use n=4 as a representative example.
    """
    n = 4

    print("The prefactor 'c_n' is derived from the standard formula for the n-th virial coefficient B_n:")
    print("B_n = -((n-1)/n!) * Integral(Sum over biconnected graphs on n vertices)")
    print("\nThe 'fully f-connected' diagram corresponds to the complete graph K_n.")
    print("Its contribution to B_n is isolated to find the prefactor c_n.")
    print("This gives the formula: c_n = -(n-1) / n!")
    print("\n--- Calculation for n = 4 ---")

    # Step 1: Calculate the numerator of the formula.
    numerator = n - 1
    print(f"The numerator is (n-1): {n} - 1 = {numerator}")

    # Step 2: Calculate the denominator of the formula.
    denominator = math.factorial(n)
    print(f"The denominator is n!: {n}! = {denominator}")

    # Step 3: Print the final equation with the calculated numbers.
    print(f"The final equation for c_4 is: c_4 = -({numerator}) / ({denominator})")

    # To provide a cleaner final answer, simplify the fraction.
    common_divisor = math.gcd(numerator, denominator)
    simplified_numerator = numerator // common_divisor
    simplified_denominator = denominator // common_divisor

    print(f"As a simplified fraction: c_4 = -{simplified_numerator}/{simplified_denominator}")
    
    # Also calculate and print the float value.
    c_n_float = -numerator / denominator
    print(f"As a decimal: c_4 = {c_n_float}")

if __name__ == "__main__":
    calculate_prefactor_cn()