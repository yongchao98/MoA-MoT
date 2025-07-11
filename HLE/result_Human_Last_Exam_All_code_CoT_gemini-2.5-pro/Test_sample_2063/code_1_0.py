# Define the base of the probability (number of choices for the remainder)
base = 8

# The number of "successful" outcomes for the remainder of a single base count.
# The remainder c_x mod 8 can be any value from 0 to 7.
# The collapse condition is c_x mod 8 = x.
# So for each x, there is 1 failing value and (8-1)=7 successful values.
successful_outcomes = base - 1

# The number of bases we need to check (from x=0 to x=7)
num_bases = 8

# The limiting probability is the product of the individual probabilities
# for each base, which is (7/8) multiplied by itself 8 times.
numerator = successful_outcomes ** num_bases
denominator = base ** num_bases

# Calculate the final probability
result = numerator / denominator

# Print the final equation with each number, as requested.
print(f"The closed form expression for the limiting probability is ({successful_outcomes}/{base})^{num_bases}.")
print(f"This is calculated as: {numerator} / {denominator}")
print(f"The final limiting probability is: {result}")
