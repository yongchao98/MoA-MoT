import numpy as np

def get_valuation(n, p):
    """Computes the p-adic valuation of an integer n."""
    if n == 0:
        return float('inf')
    val = 0
    while n % p == 0:
        val += 1
        n //= p
    return val

def poly_resultant(p, q):
    """Computes the resultant of two polynomials p and q."""
    # p, q are lists of coefficients in descending order of power
    deg_p = len(p) - 1
    deg_q = len(q) - 1
    if deg_p < 0 or deg_q < 0:
        return 0
    
    # Sylvester matrix
    sylvester_matrix = np.zeros((deg_p + deg_q, deg_p + deg_q), dtype=np.int64)
    
    for i in range(deg_q):
        sylvester_matrix[i, i:i+deg_p+1] = p
        
    for i in range(deg_p):
        sylvester_matrix[i+deg_q, i:i+deg_q+1] = q
        
    return int(round(np.linalg.det(sylvester_matrix)))

# Define the polynomial f(x) = 2*x^5 + 2*x^3 + 1
f = [2, 0, 2, 0, 0, 1]
# Define its derivative f'(x) = 10*x^4 + 6*x^2
f_prime = [10, 0, 6, 0, 0]

# Degree and leading coefficient of f(x)
n = len(f) - 1
a_n = f[0]

# Calculate the resultant
res = poly_resultant(f, f_prime)

# Calculate the discriminant
# The sign is (-1)^(n(n-1)/2) = (-1)^10 = 1
discriminant = res // a_n

# Get the 2-adic valuation of the discriminant
val_discriminant = get_valuation(discriminant, 2)

# Calculate the naive valuation based on root analysis
naive_val_discriminant = get_valuation(a_n**(2*n-2), 2) - 4

# The thickness is the difference
thickness = val_discriminant - naive_val_discriminant

print(f"The polynomial is f(x) = 2*x^5 + 2*x^3 + 1")
print(f"The discriminant of f(x) is {discriminant}")
print(f"The 2-adic valuation of the discriminant is v_2(Δ) = {val_discriminant}")
print(f"The 'naive' valuation from root analysis is {naive_val_discriminant}")
print(f"The thickness δ is the difference: δ = v_2(Δ) - v_2(Δ_naive)")
print(f"δ = {val_discriminant} - {naive_val_discriminant} = {thickness}")
