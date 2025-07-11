import math

# The problem can be solved analytically through integral calculus. The final result for
# the difference between the probability of an electron escaping through the hypotenuse
# and the probability of escaping through either of the two legs is given by the formula:
# Probability_Difference = -ln(2) / pi
# where ln(2) is the natural logarithm of 2.

# First, we define the mathematical constants needed.
natural_log_of_2 = math.log(2)
pi_value = math.pi

# Now, we calculate the final result using the formula.
result = -natural_log_of_2 / pi_value

print("The analytical solution for the difference in probabilities is given by the equation:")
print("Probability_Difference = -ln(2) / pi\n")

print("The values used in the equation are:")
print(f"ln(2) = {natural_log_of_2}")
print(f"pi = {pi_value}\n")

print("Substituting the numbers into the equation:")
print(f"Probability_Difference = -{natural_log_of_2} / {pi_value}\n")

print("The final calculated difference in probabilities is:")
print(result)