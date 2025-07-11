import sys

def solve():
    """
    Computes the value of the limit lim_{m -> infinity} ln(f(m))/ln(m).
    The user should provide the integer k >= 2 as a command-line argument.
    """
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py <k>")
        print("Please provide an integer k >= 2 as an argument.")
        return

    try:
        k = int(sys.argv[1])
        if k < 2:
            raise ValueError()
    except ValueError:
        print("Invalid input. Please provide an integer k >= 2.")
        return

    # Based on the reasoning that the worst-case matrix is a dense square block of ones,
    # the exponent in the asymptotic behavior of f(m) is 1 - 1/(2*k).
    # We will now calculate this value.

    # The components of the final formula are:
    a = 1
    b = 1
    d = 2

    print(f"For k = {k}, the formula for the limit is c = a - b / (d * k)")
    print(f"Plugging in the numbers, we get:")
    print(f"c = {a} - {b} / ({d} * {k})")

    # Perform the calculation
    result_numerator = a * d * k - b
    result_denominator = d * k

    # Simplify the fraction by finding the greatest common divisor
    def gcd(x, y):
        while y:
            x, y = y, x % y
        return abs(x)

    common_divisor = gcd(result_numerator, result_denominator)
    simplified_num = result_numerator // common_divisor
    simplified_den = result_denominator // common_divisor
    
    if simplified_den == 1:
        final_result_str = f"{simplified_num}"
    else:
        final_result_str = f"{simplified_num}/{simplified_den}"

    print(f"The result is {result_numerator}/{result_denominator}, which simplifies to {final_result_str}.")

if __name__ == "__main__":
    solve()