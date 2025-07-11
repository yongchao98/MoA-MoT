import sympy

def solve_curve_problem():
    """
    Analyzes a list of hyperelliptic curves to find which one has good reduction
    for all prime numbers greater than 2.

    This script requires the 'sympy' library. You can install it using:
    pip install sympy
    """

    print("Analyzing which curve has good ordinary reduction above 2.")
    print("A necessary condition for this is having 'good reduction' for all primes p > 2.")
    print("This means the discriminant of the polynomial f(x) in z^2 = f(x) must only have 2 as a prime factor.\n")

    x = sympy.Symbol('x')
    
    # The polynomials f(x) from the answer choices z^2 = f(x)
    # The coefficients and powers for each equation are explicitly used below.
    curves = {
        'A': x**5 + 3,
        'B': x**5 - 1,
        'C': x**6 - 1,
        'D': 2*x**5 + 2*x**3 + 1,
        'E': 4*x**5 + 4*x**3 + x**2 + 4*x
    }
    
    equations = {
        'A': "z^2 = 1*x^5 + 3",
        'B': "z^2 = 1*x^5 - 1",
        'C': "z^2 = 1*x^6 - 1",
        'D': "z^2 = 2*x^5 + 2*x^3 + 1",
        'E': "z^2 = 4*x^5 + 4*x^3 + 1*x^2 + 4*x"
    }

    correct_choice = None

    for label, f in curves.items():
        # Retrieve the full equation string for printing
        equation = equations[label]
        print(f"--- Analyzing Curve {label}: {equation} ---")

        # 1. Calculate the discriminant
        disc = sympy.discriminant(f, x)
        print(f"The discriminant is: {disc}")

        if disc == 0:
            print("Result: The polynomial has multiple roots, leading to bad reduction.")
            has_good_reduction_above_2 = False
        else:
            # 2. Find prime factors of the discriminant's absolute value
            factors = list(sympy.factorint(abs(disc)).keys())
            print(f"Prime factors of the discriminant are: {factors}")

            # 3. Check for prime factors greater than 2
            bad_primes = [p for p in factors if p > 2]
            
            if not bad_primes:
                print("Result: The discriminant has no prime factors > 2. This curve has good reduction above 2.")
                has_good_reduction_above_2 = True
                correct_choice = label
            else:
                print(f"Result: The discriminant is divisible by {bad_primes}. The curve has bad reduction at these primes.")
                has_good_reduction_above_2 = False
        
        print("-" * (len(equation) + 21) + "\n")

    if correct_choice:
        print(f"Conclusion: Curve {correct_choice} is the only option with good reduction for all primes p > 2.")
        print("Therefore, it is the answer.")
    else:
        print("Conclusion: No single curve could be identified based on the good reduction criterion.")

if __name__ == '__main__':
    solve_curve_problem()