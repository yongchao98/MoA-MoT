# The analysis of the RASPy code reveals that for the given inputs,
# the function 'v' will perform a simple addition of the two numbers
# provided in the input string. The complex "pwned" detection logic
# is not triggered.

# The task is to find the result for two separate calculations.

# First calculation, for the input '734107+4295754'.
num1_first_eq = 734107
num2_first_eq = 4295754
result1 = num1_first_eq + num2_first_eq

# Second calculation, for the input '5429141+142196'.
num1_second_eq = 5429141
num2_second_eq = 142196
result2 = num1_second_eq + num2_second_eq

# The final output is requested in the format "output1;output2".
print(f"{result1};{result2}")