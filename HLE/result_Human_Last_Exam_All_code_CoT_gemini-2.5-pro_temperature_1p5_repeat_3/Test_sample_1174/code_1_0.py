import sys

# Define the exponents based on the physical derivation
n1 = -1
n2 = 1
n3 = -3
n4 = 0
n5 = 1
n6 = 0.5

# Calculate each term in the sum
term1 = 1 * n1
term2 = 2 * n2
term3 = 3 * n3
term4 = 4 * n4
term5 = 5 * n5
term6 = 6 * n6

# Calculate the final sum
total_sum = term1 + term2 + term3 + term4 + term5 + term6

# Print the equation with the values
print(f"The calculation is:")
# Use a custom function to format positive numbers with a leading '+'
def format_num(n):
    return f"+ {n}" if n >= 0 else f"- {-n}"

# We need to output the equation exactly as requested in the prompt
# "output each number in the final equation"
# Let's print the full expression before evaluating
print(f"1*({n1}) + 2*({n2}) + 3*({n3}) + 4*({n4}) + 5*({n5}) + 6*({n6})")
print(f"= {term1} {format_num(term2)} {format_num(term3)} {format_num(term4)} {format_num(term5)} {format_num(term6)}")
print(f"= {total_sum}")

# Use this to output the final numerical answer for the platform
sys.stdout.write(f'<<<{total_sum}>>>')