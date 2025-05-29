# Calculate the expected value of an eight-sided die
outcomes = range(1, 9)  # Possible outcomes from 1 to 8
probability = 1 / 8  # Probability of each outcome

# Calculate the expected value
expected_value = sum(outcome * probability for outcome in outcomes)
print(expected_value)