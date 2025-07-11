import sympy
import math

def analyze_curve(option_letter, poly, x_symbol):
    """
    Analyzes a hyperelliptic curve z^2 = poly to determine if it has good
    reduction for all primes p > 2.
    """
    print(f"--- Analyzing Option {option_letter}: z^2 = {poly} ---")

    # 1. Find the non-singular model by removing square factors from the polynomial.
    # The non-singular model is z^2 = g(x), where g(x) is the square-free part of poly.
    try:
        sqf_list = sympy.sqf_list(poly)[1]
        square_free_part = sympy.prod(p for p, e in sqf_list if e % 2 != 0)
        
        has_square_factors = any(e >= 2 for p, e in sqf_list)
        if has_square_factors:
            print(f"The polynomial has square factors. The curve is birationally equivalent to the non-singular model: z^2 = {square_free_part}")
        else:
            print("The polynomial is square-free, so the model is already non-singular.")
    except Exception as e:
        print(f"Could not process the polynomial: {e}")
        return

    # 2. Determine the genus of the non-singular model.
    # Genus g = floor((degree - 1) / 2) for a hyperelliptic curve.
    degree = sympy.degree(square_free_part, x_symbol)
    if degree < 1:
        # Not a standard curve for this analysis
        genus = -1
        print("The non-singular model is not a standard hyperelliptic curve (degree < 1).")
        return
    else:
        genus = math.floor((degree - 1) / 2)
    print(f"The polynomial of the non-singular model has degree {degree}, so the curve has genus {genus}.")

    # 3. Check for the good reduction condition.
    if genus == 0:
        print("A genus 0 curve defined over the rational numbers has good reduction at all primes.")
        print("Conclusion: This curve HAS good ordinary reduction above 2.")
        return

    # For genus > 0, we check the discriminant.
    # The discriminant must only be divisible by 2 for good reduction at p > 2.
    # We first create a monic polynomial with integer coefficients to compute the discriminant correctly.
    integral_poly, _ = square_free_part.as_poly(x_symbol).clear_den()
    disc = sympy.discriminant(integral_poly)

    if disc == 0:
        print("The discriminant is 0. This is unexpected for a square-free polynomial.")
        return
        
    abs_disc = abs(int(disc))
    prime_factors = list(sympy.factorint(abs_disc).keys())
    odd_prime_factors = [p for p in prime_factors if p > 2]

    print(f"The discriminant of the polynomial is {disc}.")
    print(f"The prime factors of its absolute value are: {prime_factors}.")

    if not odd_prime_factors:
        # This case is not reached by any of the options with genus > 0
        print("The discriminant is a power of 2. The curve has good reduction for all primes p > 2.")
    else:
        print(f"The discriminant is divisible by the odd prime(s) {odd_prime_factors}.")
        print("Conclusion: This curve has bad reduction at primes > 2 and does NOT satisfy the condition.")


def solve_task():
    """
    Main function to analyze all curves and find the one with good ordinary reduction above 2.
    """
    x = sympy.Symbol('x')

    # The equations from the answer choices
    equations = {
        'A': x**5 + 3,
        'B': x**5 - 1,
        'C': x**6 - 1,
        'D': 2*x**5 + 2*x**3 + 1,
        'E': 4*x**5 + 4*x**3 + x**2 + 4*x
    }
    
    # Analyze each equation
    for option, poly in equations.items():
        analyze_curve(option, poly, x)
        print("\n" + "="*50 + "\n")

if __name__ == '__main__':
    solve_task()