import math

# Define the constants
sqrt5 = math.sqrt(5)
lambda1 = (1 + sqrt5) / 2
lambda2 = (1 - sqrt5) / 2
t_val = 10 * math.log(10) # This is ln(10^10)

# The expression to calculate is -(phi_0(t_val)) + 2/sqrt(5)
# Based on our analysis, this simplifies to (1/sqrt(5)) * (exp(lambda1*t) + exp(lambda2*t))
# Let's define the components of the final equation
c0 = 1
c1 = 5
c2 = 1
c3 = 5
c4 = 2
c5 = 10
c6 = 10
c7 = 1
c8 = 5
c9 = 10
c10 = 1
c11 = 5
c12 = 2

# Calculate the final value
# value = (1/sqrt5) * (math.exp(lambda1 * t_val) + math.exp(lambda2 * t_val))
# This is equivalent to: (1/sqrt5) * (10**(10*lambda1) + 10**(10*lambda2))
# Let's use the exponential form for precision with large numbers
term1 = math.exp(lambda1 * t_val)
term2 = math.exp(lambda2 * t_val)
value = (1 / sqrt5) * (term1 + term2)

# Print the final equation with its components
print(f"The final expression to calculate is of the form:")
print(f"Value = (1/sqrt({c1})) * (exp((({c2}+{c3}**0.5)/{c4}) * {c5}*ln({c6})) + exp((({c7}-{c8}**0.5)/{c12}) * {c9}*ln({c6})))")
print("\nWhich simplifies to:")
print(f"Value = (1/math.sqrt(5)) * (math.exp(lambda1 * t) + math.exp(lambda2 * t)) where t = 10*ln(10)")

# Print the final calculated value
print(f"\nCalculated Value: {value}")
