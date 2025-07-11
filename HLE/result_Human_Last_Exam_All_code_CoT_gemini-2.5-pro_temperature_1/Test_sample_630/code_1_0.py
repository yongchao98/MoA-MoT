import sympy

# Define symbols for the analysis
# T is the number of samples
# D is the radius of the feasible set
T, D = sympy.symbols('T D', positive=True)

# Define the exponents for our error model:
# Variance is proportional to D^b / T^a
# Bias is proportional to D^(-c)
a = 1  # Exponent of T in Variance term
b = 1  # Exponent of D in Variance term
c = 2  # Exponent of D in Bias term

print(f"Modeling Error = (D^{b}/T^{a}) + D^{-c}")
print(f"Variance term exponent for T: a = {a}")
print(f"Variance term exponent for D: b = {b}")
print(f"Bias term exponent for D: c = {c}")
print("-" * 20)

# To balance the terms, we set them equal (ignoring constants)
# D^b / T^a = D^(-c)
# T^a = D^(b+c)
# D = T^(a / (b+c))
# The exponent of T in the optimal D is p
p = sympy.Rational(a, b + c)

print(f"The optimal radius D scales with T as: D ~ T^p")
print(f"where p = a / (b + c) = {a} / ({b} + {c}) = {p}")
print("-" * 20)

# The optimal rate is the value of the error at this optimal D.
# We can calculate it from the bias term D^(-c)
# Rate ~ (T^p)^(-c) = T^(-c*p)
final_exponent = -c * p

print(f"The convergence rate is then O(T^q)")
print(f"where q = -c * p = -{c} * {p} = {final_exponent}")
print("-" * 20)
print(f"So, the optimal rate of convergence is Theta(1/T^{sympy.denom(final_exponent)}/{sympy.numer(final_exponent)})")
# Final answer requires outputting the numbers in the final equation
print(f"The final equation for the rate is O(1/T^(2/3))")
