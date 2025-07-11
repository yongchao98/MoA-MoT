# The task is to find the "doppelgänger" (duplicate) number in a list.
# The list is assumed to contain numbers from 1 to N, with one number repeated.
numbers = [1, 2, 3, 4, 5, 6, 7, 8, 4]

# 1. Calculate the actual sum of all numbers in the list.
actual_sum = sum(numbers)

# 2. Determine the value of 'N' for the ideal, unique sequence.
# This would be the largest number in the list if there were no gaps.
# We find the max number in the list to define the sequence 1 to N.
n = max(numbers)

# 3. Calculate the sum of the expected unique sequence from 1 to N.
# The formula for the sum of an arithmetic series is N * (N + 1) / 2.
expected_sum = n * (n + 1) // 2

# 4. The duplicate number is the difference between the actual sum and the expected sum.
duplicate_number = actual_sum - expected_sum

# 5. Print the final equation with each number.
# The logic is: (Sum of list with duplicate) - (Sum of list without duplicate) = Duplicate
print("Equation to find the doppelgänger number:")
print(f"{actual_sum} - {expected_sum} = {duplicate_number}")
