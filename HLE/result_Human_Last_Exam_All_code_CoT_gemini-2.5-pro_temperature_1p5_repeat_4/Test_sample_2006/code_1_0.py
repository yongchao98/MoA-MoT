# The first number is derived from "An August tempest".
# The month August has 31 days.
num1 = 31

# The second number comes from "Twice fifteen".
# 2 * 15 = 30.
num2 = 30

# The operation is hinted at by "The Bays" (Bytes) and "winds of change".
# This suggests a bitwise XOR operation (^).
op_symbol = "^"

# We perform the bitwise XOR calculation.
result = num1 ^ num2

# The haiku states "A divine one yields", meaning the result is 1.
# Our calculation confirms this: 31 (0b11111) XOR 30 (0b11110) = 1 (0b00001).
# Now, we print the final equation, showing each number.
print(f"{num1} {op_symbol} {num2} = {result}")