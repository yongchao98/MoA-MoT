import math

# Based on the derivation, the original integral simplifies to:
# I = integral from 0 to pi of log(sin(x)) dx - integral from 0 to pi of log(x) dx

# The value of the first integral is a known result.
# I1 = integral from 0 to pi of log(sin(x)) dx = -pi * log(2)
I1_val = -math.pi * math.log(2)

# The value of the second integral is found by integration.
# I2 = integral from 0 to pi of log(x) dx = pi*log(pi) - pi
I2_val = math.pi * math.log(math.pi) - math.pi

# The final value of the original integral is I = I1 - I2
final_value = I1_val - I2_val

# Print the final equation with numerical values
print("The value of the integral is I = I1 - I2")
print(f"I = ({I1_val}) - ({I2_val})")
print(f"The final value is: {final_value}")

# To show each number in the final equation as requested:
pi_val = math.pi
log2_val = math.log(2)
logpi_val = math.log(math.pi)
print("\nThe final equation with evaluated constants is:")
print(f"I = (-{pi_val} * {log2_val}) - ({pi_val} * {logpi_val} - {pi_val})")
