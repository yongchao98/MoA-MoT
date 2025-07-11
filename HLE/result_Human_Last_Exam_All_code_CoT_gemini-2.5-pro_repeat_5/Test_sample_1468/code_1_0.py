import sympy

# Define the symbols
N = sympy.Symbol('N', positive=True, integer=True)
d_prime = sympy.Symbol("d'", positive=True, integer=True)
q = sympy.Symbol('q', positive=True, integer=True)

# The lower bound for m is Omega(N * d' / q)
lower_bound = (N * d_prime) / q

print("The asymptotic lower bound for m is Omega(f(N, d', q)), where f(N, d', q) is:")
sympy.pprint(lower_bound)

# Example values
N_val = 1000
d_prime_val = 50
q_val = 10

# Calculate the value of the expression for the example
example_value = (N_val * d_prime_val) / q_val

print(f"\nFor example, if N = {N_val}, d' = {d_prime_val}, and q = {q_val}, the expression evaluates to:")
print(f"{N_val} * {d_prime_val} / {q_val} = {example_value}")
