# The problem asks for the maximal guaranteed probability of success for Alice.

# We determined the optimal strategy for Alice is to open 19 boxes randomly
# and use the minimum and maximum of the observed values to form her
# guessing interval.

# The total number of distinct numbers in the boxes.
total_numbers = 20

# Alice's strategy fails only if the number in the single closed box is either
# the absolute minimum of all 20 numbers or the absolute maximum.
# Let the sorted numbers be y_1, y_2, ..., y_20.
# Failure occurs if the closed box contains y_1 or y_20.
losing_outcomes = 2

# The strategy succeeds if the closed box contains any of the other numbers,
# from y_2 to y_19.
winning_outcomes = total_numbers - losing_outcomes

# The probability of success is the ratio of winning outcomes to the total
# number of possible outcomes for the value in the closed box. Alice's
# random choice of which box to leave closed makes each outcome equally likely.

# The equation for the probability p is:
# p = winning_outcomes / total_numbers
print("The calculation for the maximal probability p is as follows:")
print(f"Total number of items in the sequence = {total_numbers}")
print(f"Number of 'losing' items (the global min and max) = {losing_outcomes}")
print(f"Number of 'winning' items = {total_numbers} - {losing_outcomes} = {winning_outcomes}")
print(f"The equation for the probability is: p = {winning_outcomes} / {total_numbers}")

# Calculate the final value
final_probability = winning_outcomes / total_numbers

print(f"The maximal probability p is {final_probability}, which is 9/10.")
