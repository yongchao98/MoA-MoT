import math

# Step 1: Define the values needed for the calculation
# F3(1) is ln(1 + e^(sin(1)))
# F3(0) is ln(1 + 1) = ln(2)

# Step 2: Calculate F3(1)
f3_at_1 = math.log(1 + math.exp(math.sin(1)))

# Step 3: Calculate F3(0)
f3_at_0 = math.log(2)

# Step 4: Calculate the value of the integral V
# V = ln(F3(1)) - ln(F3(0)) = ln(F3(1) / F3(0))
V = math.log(f3_at_1 / f3_at_0)

# Step 5: Calculate the final requested value
final_answer = round(10000 * V)

# Step 6: Print the output as requested
# "output each number in the final equation!"
# The final equation is ln(F3(1)) - ln(F3(0))
print(f"The evaluation of the integral results in the equation: ln(F3(1)) - ln(F3(0))")
print(f"The number for F3(1) is: {f3_at_1}")
print(f"The number for F3(0) is: {f3_at_0}")
print(f"The value of the integral V is: {V}")
print(f"The closest integer to 10000*V is: {final_answer}")