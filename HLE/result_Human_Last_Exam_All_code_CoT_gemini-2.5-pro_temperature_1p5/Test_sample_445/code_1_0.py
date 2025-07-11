# The total number of boxes.
n = 20

# Alice's strategy is to open n-1 boxes.
# This means there is 1 box left closed.
# The number in the closed box could be any of the n numbers
# with equal probability, based on Alice's random choice of which box to close.

# Let the n distinct numbers be sorted: y_1 < y_2 < ... < y_n.

# Alice wins if the number in the closed box is NOT the absolute minimum (y_1)
# and NOT the absolute maximum (y_n).
# This is because if the number is y_j where 1 < j < n, the numbers she observes
# will include y_1 and y_n. Her interval (V_min, V_max) will be (y_1, y_n),
# which contains y_j.
# If the number is y_1 or y_n, it will fall outside the observed range.

# Total number of possible values for the closed box.
total_cases = n

# Number of values for the closed box that lead to a loss.
# This happens if the box contains the minimum (y_1) or maximum (y_n) value.
losing_cases = 2

# Number of values for the closed box that lead to a win.
winning_cases = total_cases - losing_cases

# The probability of success is the ratio of winning cases to total cases.
# We use integer arithmetic for the numerator and denominator first.
p_numerator = winning_cases
p_denominator = total_cases

# Calculate the decimal value of the probability.
p_value = p_numerator / p_denominator

# Print the explanation and the final equation.
print("The maximal probability p can be found using the optimal strategy.")
print("Optimal strategy: Alice opens 19 of the 20 boxes.")
print("She wins if the number in the single closed box is not the minimum or maximum of all 20 numbers.")
print(f"Total numbers: n = {n}")
print(f"Number of 'extreme' values (min/max) that cause a loss: 2")
print(f"Number of 'inner' values that cause a win: n - 2 = {n} - 2 = {winning_cases}")
print("\nThe probability of success p is the ratio of winning cases to total cases.")
print("p = (n - 2) / n")
print(f"p = ({n} - 2) / {n} = {p_numerator}/{p_denominator}")
print(f"p = {p_value}")