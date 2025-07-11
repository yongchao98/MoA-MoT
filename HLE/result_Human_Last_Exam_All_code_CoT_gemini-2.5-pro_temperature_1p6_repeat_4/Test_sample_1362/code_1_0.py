import sympy

def analyze_curve(label, f_str):
    """
    Analyzes the properties of the curve z^2 = f(x).
    
    Args:
        label (str): The label of the curve (e.g., 'A').
        f_str (str): The string representation of the polynomial f(x).
    """
    x = sympy.Symbol('x')
    f = sympy.sympify(f_str)
    
    print(f"--- Analyzing Curve {label}: z^2 = {f} ---")
    
    # Calculate the discriminant
    try:
        disc = sympy.discriminant(f)
        print(f"Discriminant: {disc}")
        
        # Get the prime factorization of the absolute value of the discriminant
        if disc != 0:
            prime_factors = sympy.factorint(abs(disc))
            print(f"Prime factors of discriminant: {prime_factors}")
            
            odd_prime_factors = {p for p in prime_factors if p > 2}
            if not odd_prime_factors:
                print("Result: This curve HAS good reduction for all primes p > 2.")
            else:
                print(f"Result: This curve has BAD reduction at odd prime(s) {sorted(list(odd_prime_factors))}.")
        else:
            print("Result: Discriminant is 0, indicating repeated roots and bad reduction.")

    except Exception as e:
        print(f"Could not compute discriminant: {e}")
    
    # Check for reducibility over the rational numbers
    try:
        factors = sympy.factor(f)
        if factors == f:
            print("The polynomial f(x) is irreducible over the rationals.")
        else:
            print(f"The polynomial f(x) is reducible over the rationals: {factors}")
    except Exception as e:
        print(f"Could not factor the polynomial: {e}")
    print("-" * (len(label) + 26))


def main():
    curves = {
        'A': 'x**5 + 3',
        'B': 'x**5 - 1',
        'C': 'x**6 - 1',
        'D': '2*x**5 + 2*x**3 + 1',
        'E': '4*x**5 + 4*x**3 + x**2 + 4*x'
    }
    
    for label, f_str in curves.items():
        analyze_curve(label, f_str)

if __name__ == "__main__":
    main()