import math

def solve():
    """
    Calculates the value of P(n) for a user-specified integer n.
    """
    try:
        n_str = input("Enter a positive integer value for n: ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    # Calculate components of the formula
    base1 = 2 * math.pi
    c1 = 2
    
    exponent1_num = n * (n - 1)
    exponent1_den = 4
    exponent1 = exponent1_num / exponent1_den
    
    base2 = math.factorial(n)
    
    exponent2_num = -1
    exponent2_den = 2
    exponent2 = exponent2_num / exponent2_den

    # Calculate the final result
    result = (base1 ** exponent1) * (base2 ** exponent2)

    # Output the formula and the result
    print("\nThe closed-form formula is: P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print(f"\nFor n = {n}:")
    print(f"The first term in the numerator is ({c1} * pi)^({exponent1_num} / {exponent1_den})")
    print(f"The term in the denominator is sqrt({base2})")
    print(f"\nP({n}) = ({base1:.4f})^({exponent1}) / {math.sqrt(base2):.4f}")
    print(f"P({n}) = {result}")

solve()
