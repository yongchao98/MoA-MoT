# The total number of boxes in the sequence.
N = 20

# Alice's strategy involves opening N-1 boxes and guessing the Nth.
# Her guess wins unless the unopened box contains the global minimum or the global maximum
# of all N numbers.
# Therefore, there are 2 losing cases based on the rank of the number in the unopened box.
num_losing_cases = 2

# The total number of possibilities for the rank of the number in the unopened box is N.
num_total_cases = N

# The number of winning cases is the total number of cases minus the losing cases.
num_winning_cases = num_total_cases - num_losing_cases

# The probability of success is the ratio of winning cases to total cases.
# We will use integer division to get the simplified fraction's numerator and denominator.
import math
common_divisor = math.gcd(num_winning_cases, num_total_cases)
numerator = num_winning_cases // common_divisor
denominator = num_total_cases // common_divisor

print(f"Total number of boxes (N): {num_total_cases}")
print(f"Number of losing cases (if the unopened box contains the overall min or max): {num_losing_cases}")
print(f"Number of winning cases = {num_total_cases} - {num_losing_cases} = {num_winning_cases}")
print(f"The maximal probability p = Winning Cases / Total Cases")
print(f"p = {num_winning_cases} / {num_total_cases}")
print(f"p = {numerator} / {denominator}")
