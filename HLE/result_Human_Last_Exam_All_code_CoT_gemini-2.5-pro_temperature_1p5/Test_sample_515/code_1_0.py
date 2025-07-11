import math

# Constants
gamma = 0.5772156649
log2 = math.log(2)
log3600 = math.log(3600)

# Potential kernel at x0 = (0,1) is a(1,0)
# a(1,0) = 2/pi
a_x0 = 2 / math.pi

# Asymptotic formula for potential kernel at z = (3600,0)
# a(z) ~ (2/pi) * ln(||z||) + (2*gamma + ln(8)) / pi
# The constant term is (2*gamma + 3*ln(2))/pi
a_z = (2 / math.pi) * log3600 + (2 * gamma + 3 * log2) / math.pi

# The probability p is the ratio a(x0) / a(z)
# The 2/pi term cancels out from the numerator and the first part of the denominator.
# p ~ 1 / (ln(3600) + gamma + 1.5*ln(2))
denominator = log3600 + gamma + 1.5 * log2
prob = 1 / denominator

print(f"The calculation is based on the formula P ~= a(x_0)/a(z).")
print(f"a(x_0) = a(0,1) = a(1,0), which is exactly 2/pi.")
print(f"a(z) = a(3600,0) is approximated by the asymptotic formula (2/pi)*ln(3600) + (2*gamma + ln(8))/pi.")
print(f"The probability is approximately 1 / (ln(3600) + gamma + 1.5*ln(2)).")

term1 = log3600
term2 = gamma
term3 = 1.5 * log2

# Using print to show the final equation with computed values
# Using f-string formatting to control the number of digits for clarity
print(f"p ~= 1 / ({term1:.4f} + {term2:.4f} + {term3:.4f})")
print(f"p ~= 1 / {denominator:.4f}")
print(f"p ~= {prob:.4f}")
# The result should be given with two significant digits
print(f"The approximate probability is {prob:.2g}.")
print("Final Answer in numerical form for evaluation:")
print(f"{prob:.2g}")