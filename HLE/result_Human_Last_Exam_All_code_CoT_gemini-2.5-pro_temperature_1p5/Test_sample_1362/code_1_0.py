# This script requires the sympy library.
# You can install it by running: pip install sympy

import sympy

def find_curve_with_good_reduction():
    """
    Analyzes a list of hyperelliptic curves to determine which one has 
    good reduction for all prime numbers p > 2.
    """
    # Define the symbol 'x' for our polynomials
    x = sympy.Symbol('x')

    # Store the polynomials f(x) from each option z^2 = f(x) in a dictionary
    # along with a string representation of the equation.
    curves = {
        'A': (x**5 + 3, "z^2 = x**5 + 3"),
        'B': (x**5 - 1, "z^2 = x**5 - 1"),
        'C': (x**6 - 1, "z^2 = x^6 - 1"),
        'D': (2*x**5 + 2*x**3 + 1, "z^2 = 2*x**5 + 2*x**3 + 1"),
        'E': (4*x**5 + 4*x**3 + x**2 + 4*x, "z^2 = 4*x**5 + 4*x**3 + x**2 + 4*x")
    }

    print("Analyzing which curve has good reduction for all primes p > 2.")
    print("This requires that the discriminant of the polynomial f(x) in z^2=f(x) is not divisible by any prime other than 2.\n")

    correct_option = None

    for option, (poly, eqn_str) in curves.items():
        print(f"--- Checking Curve {option}: {eqn_str} ---")
        
        # 1. Calculate the discriminant of the polynomial
        # The int() conversion handles sympy's integer type
        disc = int(sympy.discriminant(poly, x))
        print(f"The discriminant is: {disc}")
        
        # 2. Find the prime factors of the discriminant's absolute value
        if disc == 0:
            prime_factors = []
            print("Discriminant is 0, which means the polynomial has repeated roots.")
        else:
            factors_dict = sympy.factorint(abs(disc))
            prime_factors = sorted(factors_dict.keys())
            print(f"The prime factors of the discriminant are: {prime_factors}")

        # 3. Check if any prime factor is greater than 2
        primes_above_2 = [p for p in prime_factors if p > 2]
        
        if not primes_above_2:
            print("Result: This curve has GOOD reduction for all primes p > 2.")
            if correct_option is None: # Store the first one found
                 correct_option = option
        else:
            print(f"Result: This curve has BAD reduction at prime(s) {primes_above_2}.")
        
        print("-" * 50)
        print()

    if correct_option:
        print("\n--- CONCLUSION ---")
        print(f"Curve {correct_option} is the correct answer.")
        final_eqn = curves[correct_option][1]
        print(f"It is the only choice where the discriminant's prime factors do not include any prime greater than 2.")
        print(f"The equation for the curve is: {final_eqn}")

find_curve_with_good_reduction()