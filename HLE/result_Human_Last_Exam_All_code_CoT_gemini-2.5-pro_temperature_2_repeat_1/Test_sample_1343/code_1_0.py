# --- Calculation for the first input string ---
# Extract the two numbers from the input string.
# The original problem defines sop as the full string.
sop1 = "734107+4295754"
num1_str, num2_str = sop1.split('+')

# Convert string numbers to integers for calculation.
num1 = int(num1_str)
num2 = int(num2_str)

# Perform the addition.
result1 = num1 + num2


# --- Calculation for the second input string ---
# Extract the two numbers from the input string.
sop2 = "5429141+142196"
num3_str, num4_str = sop2.split('+')

# Convert string numbers to integers for calculation.
num3 = int(num3_str)
num4 = int(num4_str)

# Perform the addition.
result2 = num3 + num4

# Print the final results in the required "output1;output2" format.
# We also print the full equations to show how the numbers were derived,
# as hinted by the instructions.
print(f"Calculation 1: {num1_str} + {num2_str} = {result1}")
print(f"Calculation 2: {num3_str} + {num4_str} = {result2}")
print(f"Final combined output: {result1};{result2}")