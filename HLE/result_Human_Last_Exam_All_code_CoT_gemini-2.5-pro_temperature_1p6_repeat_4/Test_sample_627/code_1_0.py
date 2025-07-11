# The three-twist knot is the 5_2 knot in Alexander-Briggs notation.

# The Jones polynomial V(t) for the 5_2 knot is:
# V(t) = t^-2 - t^-3 + t^-4 - t^-5 + t^-6
# We represent the polynomial by the list of its exponents.
exponents = [-2, -3, -4, -5, -6]

# Find the highest and lowest powers of t in the polynomial.
max_power = max(exponents)
min_power = min(exponents)

# Calculate the span of the polynomial.
span_t = max_power - min_power

# Apply the Morton-Franks-Williams inequality to find a lower bound for the braid index b(K).
# b(K) >= (span_t(V(K)) + 2) / 2
lower_bound = (span_t + 2) / 2

print("Step 1: Find a lower bound for the braid index.")
print(f"The Jones polynomial for the three-twist knot has exponents: {exponents}")
print(f"The maximum power is {max_power}.")
print(f"The minimum power is {min_power}.")
print(f"The span of the polynomial is max_power - min_power = {max_power} - ({min_power}) = {span_t}.")
print(f"The lower bound for the braid index is (span + 2) / 2 = ({span_t} + 2) / 2 = {int(lower_bound)}.")
print("-" * 20)

print("Step 2: Find an upper bound for the braid index.")
print("The three-twist knot can be represented as the closure of a 3-strand braid.")
print("This means its braid index must be less than or equal to 3.")
upper_bound = 3
print(f"So, braid_index <= {upper_bound}.")
print("-" * 20)

print("Step 3: Combine the bounds.")
print(f"We have found that {int(lower_bound)} <= braid_index <= {upper_bound}.")
print("Since the lower and upper bounds are both 3, the braid index must be exactly 3.")
braid_index = 3
print(f"The braid index of the three-twist knot is {braid_index}.")
print("This is the least upper bound.")
