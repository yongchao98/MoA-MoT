import sympy

def find_curve_with_good_reduction_above_2():
    """
    Analyzes hyperelliptic curves to find which has good reduction for all primes p > 2.
    
    A curve z^2 = P(x) has good reduction at a prime p if p does not divide
    the discriminant of the polynomial P(x). For the curve to have good reduction
    for all primes p > 2, the discriminant of P(x) must not be divisible by any
    odd prime. Its absolute value must therefore be a power of 2.
    
    This script calculates the discriminant for each polynomial from the options
    and checks its prime factorization to find the correct curve.
    """
    print("Analyzing curves for good reduction for all primes p > 2...")
    print("="*60)

    # Dictionary mapping option labels to their polynomials P(x)
    polynomials = {
        'A': 'x**5 + 3',
        'B': 'x**5 - 1',
        'C': 'x**6 - 1',
        'D': '2*x**5 + 2*x**3 + 1',
        'E': '4*x**5 + 4*x**3 + x**2 + 4*x'
    }
    
    x = sympy.Symbol('x')
    correct_answer = None

    for label, poly_str in polynomials.items():
        print(f"Option {label}: z^2 = {poly_str}")
        
        # Create the polynomial object from the string
        P = sympy.sympify(poly_str, locals={'x': x})
        
        # Calculate the discriminant
        disc = sympy.discriminant(P, x)
        
        # Find prime factors of the absolute value of the discriminant
        factors = sympy.factorint(abs(int(disc)))
        
        # Check if 2 is the only prime factor
        is_good_above_2 = all(p == 2 for p in factors.keys())
        
        # Output the analysis for this curve, including the required numbers
        print(f"  - The polynomial is P(x) = {P}.")
        print(f"  - The discriminant of P(x) is {disc}.")
        
        # Format the prime factorization for clear output
        factor_str_parts = [f"{p}^{e}" for p, e in sorted(factors.items())]
        factor_str = " * ".join(factor_str_parts)
        if int(disc) < 0:
            factor_str = "-1 * " + factor_str
        print(f"  - The factored discriminant is: {factor_str}")
        
        if is_good_above_2:
            print("  => Conclusion: Good reduction for all p > 2 (discriminant is a power of 2).")
            correct_answer = label
        else:
            odd_primes = [p for p in factors.keys() if p != 2]
            print(f"  => Conclusion: Bad reduction at odd prime(s) {odd_primes}.")
        print("-" * 60)

    if correct_answer:
        print(f"\nFinal Answer: The correct curve is {correct_answer}.")

# Execute the analysis
find_curve_with_good_reduction_above_2()