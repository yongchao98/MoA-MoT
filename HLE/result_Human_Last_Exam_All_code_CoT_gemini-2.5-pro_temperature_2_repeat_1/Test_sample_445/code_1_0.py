# The total number of boxes.
N = 20

# Alice's strategy is to pick one box at random to make a guess on, and open the rest.
# She wins if the number in her chosen box is not the overall minimum or maximum.

# The number of outcomes where Alice's strategy fails. This occurs if she
# happens to choose the box containing the smallest (minimum) value or the
# box containing the largest (maximum) value.
num_unfavorable_outcomes = 2

# The total number of boxes she could choose as her target.
total_outcomes = N

# The number of outcomes where her strategy succeeds.
num_favorable_outcomes = total_outcomes - num_unfavorable_outcomes

# The probability of success is the ratio of favorable outcomes to the total number of outcomes.
# We are asked to output the numbers in the final equation.

print("This script calculates the maximal probability of success for Alice.")
print("The optimal strategy guarantees a win if the single unopened box does not contain an 'extreme' value (the overall minimum or maximum).")
print("")

print(f"Total number of boxes: {total_outcomes}")
print(f"Number of boxes with extreme values: {num_unfavorable_outcomes}")
print(f"Number of boxes with non-extreme values: {num_favorable_outcomes}")
print("")
print("The final equation for the probability of success (p) is:")

# We present the final equation with the numbers, as requested.
# Equation: p = (Total Boxes - Extreme Boxes) / Total Boxes
# The parts of the equation are:
eq_numerator = num_favorable_outcomes
eq_denominator = total_outcomes
# The calculation can also be seen as 1 - P(failure)
eq_1 = 1
eq_fail_num = num_unfavorable_outcomes
eq_fail_den = total_outcomes

print(f"{eq_1} - {eq_fail_num}/{eq_fail_den} = {eq_numerator}/{eq_denominator}")

# Calculate the final decimal value
p_success = num_favorable_outcomes / total_outcomes
print(f"\nThe maximal probability p is {p_success}, which is equivalent to 9/10.")
