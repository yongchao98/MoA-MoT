import sympy
from sympy import Poly

def find_curve_with_good_ordinary_reduction():
    """
    Analyzes a list of hyperelliptic curves to find one with good ordinary
    reduction for a prime p > 2.
    """
    x = sympy.symbols('x')
    
    # The polynomials f(x) for each curve z^2 = f(x)
    curves = {
        'A': 'x**5 + 3',
        'B': 'x**5 - 1',
        'C': 'x**6 - 1',
        'D': '2*x**5 + 2*x**3 + 1',
        'E': '4*x**5 + 4*x**3 + x**2 + 4*x'
    }

    # Small primes > 2 to test for reduction
    primes_to_test = [3, 5, 7]
    winner = None

    for label, f_str in curves.items():
        f = sympy.sympify(f_str)
        f_poly = Poly(f, x)
        
        # Calculate the discriminant of the polynomial f(x)
        disc = sympy.discriminant(f_poly)

        print(f"--- Analyzing Curve {label}: z^2 = {f_str} ---")
        print(f"The discriminant is {disc}.")

        # Find the first prime p > 2 where the curve has good reduction
        for p in primes_to_test:
            # 1. Good Reduction Check: Is the discriminant divisible by p?
            if disc % p == 0:
                print(f"At p={p}, the curve has BAD reduction because {p} divides the discriminant.")
                continue
            
            print(f"At p={p}, the curve has GOOD reduction.")
            
            # 2. Ordinary Reduction Check
            exponent = (p - 1) // 2
            target_power = p - 1
            
            # To be clear, we calculate the full integer coefficient first
            full_powered_poly = f_poly ** exponent
            full_coeff = full_powered_poly.coeff_monomial(x**target_power)
            
            # The actual check is done on the coefficient modulo p
            coeff_mod_p = full_coeff % p

            print(f"To check for ordinary reduction at p={p}, we test the coefficient of x^{target_power} in (f(x))^{exponent}.")
            
            # The final equation and numbers for the check
            print(f"The integer coefficient is {full_coeff}. We check this value modulo {p}.")
            print(f"The equation for the test is: {full_coeff} mod {p} = {coeff_mod_p}")

            if coeff_mod_p == 0:
                print(f"Result: The reduction is SUPERSINGULAR because the result is 0.")
            else:
                print(f"Result: The reduction is ORDINARY because the result is not 0.")
                print(f"\nCONCLUSION: Curve {label} is the answer. It has good ordinary reduction at p={p}.")
                winner = label
            
            # We only need to check the first valid prime for each curve
            break
        
        # If we found the answer, stop searching
        if winner:
            break
        
        print("-" * 55)

    if not winner:
        print("\nNone of the curves were found to have good ordinary reduction for the tested primes.")

find_curve_with_good_ordinary_reduction()