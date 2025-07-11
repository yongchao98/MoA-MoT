from fractions import Fraction

# Problem parameters
m = 4  # number of bits in B
k = 2  # number of hash functions
n = 2  # number of inputs

# Total number of hash computations
nk = n * k

# The final probability is given by the formula:
# P(B'[i]=1) = 2 * [ P(bit_j=0) - P(bit_j=0, bit_l=0) ]
# P(B'[i]=1) = 2 * [ ((m-1)/m)^(n*k) - ((m-2)/m)^(n*k) ]

# Calculate each term using the Fraction class for precision
term1_base_num = m - 1
term1_base_den = m
term1_base = Fraction(term1_base_num, term1_base_den)
term1 = term1_base ** nk

term2_base_num = m - 2
term2_base_den = m
term2_base = Fraction(term2_base_num, term2_base_den)
term2 = term2_base ** nk

# Calculate the final probability
prob = 2 * (term1 - term2)

# Print the final equation with all the numbers
print(f"The probability is calculated as follows:")
print(f"P = 2 * [ (({term1_base_num}/{term1_base_den})^{nk}) - (({term2_base_num}/{term2_base_den})^{nk}) ]")
print(f"P = 2 * [ ({term1_base})^{nk} - ({term2_base})^{nk} ]")
print(f"P = 2 * [ {term1} - {term2} ]")
print(f"P = 2 * [ {term1 - term2} ]")
print(f"P = {prob}")

# Print the final answer in the required format
print(f"\nThe final probability is {prob.numerator}/{prob.denominator}.")