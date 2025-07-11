import sys

def calculate_norm(n):
    """
    Calculates the 1-norm of the correlation matrix T for a given even integer n.
    """
    # Validate that n is a non-negative even integer
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print("Error: n must be a non-negative even integer.", file=sys.stderr)
        sys.exit(1)

    # Use Python's arbitrary-precision integers for calculations
    pow3_n = 3**n
    pow2_np1 = 2**(n + 1)

    # The simplified formula for the 1-norm of T is:
    # ||T||_1 = ((3^n - 1) * (2^(n+1) + 3) + 2) / (1 + 3^n)

    print(f"For n = {n}:")
    print(f"The 1-norm is calculated using the formula:")
    print(f"||T||_1 = ((3^n - 1) * (2^(n+1) + 3) + 2) / (1 + 3^n)")
    print(f"\nStep 1: Substitute n = {n} into the formula.")
    print(f"||T||_1 = ((3^{n} - 1) * (2^({n}+1) + 3) + 2) / (1 + 3^{n})")
    print(f"||T||_1 = (({pow3_n} - 1) * ({pow2_np1} + 3) + 2) / (1 + {pow3_n})")

    # Breaking down the calculation
    term1_in_num = pow3_n - 1
    term2_in_num = pow2_np1 + 3
    denominator = 1 + pow3_n
    print(f"\nStep 2: Evaluate the terms in the expression.")
    print(f"||T||_1 = ({term1_in_num} * {term2_in_num} + 2) / {denominator}")

    numerator_part1 = term1_in_num * term2_in_num
    print(f"\nStep 3: Perform the multiplication in the numerator.")
    print(f"||T||_1 = ({numerator_part1} + 2) / {denominator}")

    numerator = numerator_part1 + 2
    print(f"\nStep 4: Perform the addition in the numerator.")
    print(f"||T||_1 = {numerator} / {denominator}")

    # The final result should be an integer
    result = numerator // denominator
    print(f"\nStep 5: Perform the final division.")
    print(f"Final Result: {result}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <even_integer_n>", file=sys.stderr)
        sys.exit(1)
    try:
        n_val = int(sys.argv[1])
        calculate_norm(n_val)
    except ValueError:
        print("Error: Please provide a valid integer for n.", file=sys.stderr)
        sys.exit(1)
