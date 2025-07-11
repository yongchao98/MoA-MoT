import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def term_T(a, b):
    """Computes the simplified term T(a/b) from the second product."""
    x = a / b
    
    # Numerator: Gamma(x) * Gamma(x + 1/2)
    # Denominator: Gamma(2x - 1)^2
    
    # Handling negative arguments for Gamma function, as they can be poles.
    # The argument 2x-1 = (2a-b)/b is in (-1, 0) since 0 < a <= b/2.
    # Gamma function has poles at 0, -1, -2, ...
    # Since gcd(a,b)=1 and 4|b, a is odd.
    # If a=b/2, 2x-1=0, pole. But gcd(b/2, b)=b/2, so b=2, a=1. But 4|b is required.
    # So 2x-1 is never a non-positive integer.
    
    num = math.gamma(x) * math.gamma(x + 0.5)
    den = math.gamma(2 * x - 1)**2
    
    return num / den

def product_A(b):
    """Computes the product for b where 4 does not divide b."""
    prod = 1.0
    for a in range(1, b):
        if gcd(a, b) == 1:
            prod *= math.gamma(a / b)
    return prod

def product_B(b):
    """Computes the product for b where 4 divides b."""
    prod = 1.0
    for a in range(1, b // 2 + 1):
        if gcd(a, b) == 1:
            prod *= term_T(a, b)
    return prod

def calculate_P(n):
    """Calculates the value of P(n)."""
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    
    log_P = 0.0
    
    for b in range(1, n + 1):
        floor_val = n // b
        if floor_val == 0:
            continue
            
        if b % 4 != 0:
            # First product case
            log_P += floor_val * math.log(product_A(b))
        else:
            # Second product case
            log_P += floor_val * math.log(product_B(b))
            
    return math.exp(log_P)

def main():
    """Main function to get user input and print the result."""
    try:
        n_str = input("Enter a positive integer n: ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer.")
            return

        p_n = calculate_P(n)

        # The instruction "output each number in the final equation" is ambiguous.
        # For a general formula, this is complex.
        # We will output the calculated value of P(n).
        
        print(f"For n = {n}, the calculated value of P(n) is:")
        print(p_n)

    except (ValueError, TypeError):
        print("Invalid input. Please enter a positive integer.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
