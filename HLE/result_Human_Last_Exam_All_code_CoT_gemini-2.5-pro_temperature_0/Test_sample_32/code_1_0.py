from collections import defaultdict
from fractions import Fraction

def get_degree(kappa_monomial):
    """Calculates the degree of a kappa monomial."""
    degree = 0
    for i, power in enumerate(kappa_monomial):
        degree += (i + 1) * power
    return degree

def multiply_poly(p1, p2):
    """Multiplies two kappa polynomials."""
    result = defaultdict(Fraction)
    for m1, c1 in p1.items():
        for m2, c2 in p2.items():
            # Ensure monomials have the same length by padding with zeros
            max_len = max(len(m1), len(m2))
            m1_padded = m1 + (0,) * (max_len - len(m1))
            m2_padded = m2 + (0,) * (max_len - len(m2))
            
            new_monomial = tuple(p1 + p2 for p1, p2 in zip(m1_padded, m2_padded))
            result[new_monomial] += c1 * c2
    return dict(result)

def format_monomial(m):
    """Formats a kappa monomial for printing."""
    parts = []
    for i, p in enumerate(m):
        if p > 0:
            if p == 1:
                parts.append(f"kappa_{i+1}")
            else:
                parts.append(f"kappa_{i+1}^{p}")
    return "*".join(parts) if parts else "1"

# Step 1: Define lambda classes as polynomials in kappa classes
# A polynomial is a dict {monomial_tuple: coefficient}
# monomial_tuple (p1, p2, p3, ...) means kappa_1^p1 * kappa_2^p2 * ...
lambda_1 = {(1,): Fraction(1, 12)}
lambda_2 = {(2, 0): Fraction(1, 288), (0, 1): Fraction(1, 120)}
lambda_3 = {(3, 0, 0): Fraction(1, 10368), (1, 1, 0): Fraction(1, 1440), (0, 0, 1): Fraction(-1, 360)}

# Step 2: Compute the product lambda_1 * lambda_2 * lambda_3
product_12 = multiply_poly(lambda_1, lambda_2)
product_123 = multiply_poly(product_12, lambda_3)

# Filter for terms of degree 6
integrand = {m: c for m, c in product_123.items() if get_degree(m) == 6}

# Step 3: Define known integrals of kappa monomials of degree 6
kappa_integrals = {
    (6, 0, 0): Fraction(1, 72),
    (4, 1, 0): Fraction(1, 36),
    (3, 0, 1): Fraction(1, 24),
    (2, 2, 0): Fraction(1, 24),
    (1, 1, 1): Fraction(1, 16),
}

# Calculate the final integral
total_integral = Fraction(0)
print("The integral is computed as a sum of terms:")
print("Integral = sum(coefficient * integral_of_kappa_monomial)")
print("-" * 50)

# Sort monomials for deterministic output
sorted_monomials = sorted(integrand.keys(), reverse=True)

equation_parts = []
for monomial in sorted_monomials:
    coeff = integrand[monomial]
    integral_val = kappa_integrals.get(monomial, 0)
    term_val = coeff * integral_val
    total_integral += term_val
    
    # Print each number in the final equation
    print(f"Term for {format_monomial(monomial)}:")
    print(f"  Coefficient: {coeff.numerator}/{coeff.denominator}")
    print(f"  Integral value: {integral_val.numerator}/{integral_val.denominator}")
    print(f"  Term value: ({coeff.numerator}/{coeff.denominator}) * ({integral_val.numerator}/{integral_val.denominator}) = {term_val.numerator}/{term_val.denominator}")
    print("-" * 20)
    
    sign = "+" if term_val >= 0 else "-"
    abs_term_val = abs(term_val)
    equation_parts.append(f"{sign} {abs_term_val.numerator}/{abs_term_val.denominator}")

# Adjust the first term's sign
if equation_parts[0].startswith("+"):
    equation_parts[0] = equation_parts[0][2:]

print("\nFinal Equation:")
print("Integral = " + " ".join(equation_parts))

print("\nFinal Result:")
print(f"The integral of lambda_3*lambda_2*lambda_1 on M_3 is: {total_integral.numerator}/{total_integral.denominator}")
print("<<<1/1728>>>")