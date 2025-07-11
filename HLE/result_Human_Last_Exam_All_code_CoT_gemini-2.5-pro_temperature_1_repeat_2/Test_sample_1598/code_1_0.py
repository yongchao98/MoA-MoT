# The problem reduces to finding the 100th term of a minimal superincreasing sequence
# starting with 1, which follows the pattern of powers of 2.
# The sequence is 2^0, 2^1, 2^2, ...
# For n=100 edges, we need 100 distinct labels. The largest label will be 2^(100-1).

base = 2
exponent = 99

# Calculate the result
global_labeling_number = base ** exponent

# The final equation is: global_labeling_number = base ** exponent
# We print all the numbers in this equation as requested.
print(f"The calculation for the global labeling number is:")
print(f"Base: {base}")
print(f"Exponent: {exponent}")
print(f"Result: {global_labeling_number}")
