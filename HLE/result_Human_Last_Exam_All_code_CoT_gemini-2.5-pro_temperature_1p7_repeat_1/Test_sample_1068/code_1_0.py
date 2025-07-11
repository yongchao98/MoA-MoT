import math

# Set parameters based on the problem statement.
# We need to calculate l(p) for p=14.
p = 14

# The integral simplifies to a form f(a,b) where a=p and b=p+1.
a = p
b = p + 1

# The problem reduces to calculating a value from a known analytical formula
# involving the Gamma function.
# Formula: 2 * log( (Gamma(N1) * Gamma(N2)) / (Gamma(D1) * Gamma(D2)) )

# The numbers (arguments to the Gamma function) in the final equation are
# determined by a and b.
N1 = 0.5
N2 = (a + b + 1) / 2
D1 = (a + 1) / 2
D2 = (b + 1) / 2

print("The final value is calculated from the equation:")
print("Value = 2 * log( (Gamma(N1) * Gamma(N2)) / (Gamma(D1) * Gamma(D2)) )")
print("\nWhere the specific numbers for the equation are:")
print(f"N1 = {N1}")
print(f"N2 = {N2}")
print(f"D1 = {D1}")
print(f"D2 = {D2}")

# We use the math.lgamma function, which computes the natural logarithm 
# of the Gamma function, for a stable and accurate calculation.
# log( (Gamma(N1)*Gamma(N2))/(Gamma(D1)*Gamma(D2)) ) becomes
# lgamma(N1) + lgamma(N2) - lgamma(D1) - lgamma(D2).

log_gamma_part = (math.lgamma(N1) + math.lgamma(N2)) - (math.lgamma(D1) + math.lgamma(D2))

result = 2 * log_gamma_part

# The expression also has a final analytical simplification to 28 * log(2),
# which we can use to verify the result.
# We print the calculated numerical value.
print(f"\nFinal calculated result:")
print(result)