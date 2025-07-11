# The total number of chairs in the row.
N = 20

# We found that the problem reduces to finding the maximum number of non-adjacent
# people that can sit in a row.
# The formula for this is k <= (N + 1) / 2, where k is the number of people.
# Since k must be an integer, we use integer division which is equivalent to floor().
max_occupied = (N + 1) // 2

# The problem asks to output the numbers in the final equation.
print("The problem is to find the maximum number of chairs that can be occupied.")
print("This can be modeled by the formula: k = floor((N + 1) / 2)")
print(f"Given N = {N} chairs:")
print(f"k = ({N} + 1) // 2")
print(f"k = {N + 1} // 2")
print(f"k = {max_occupied}")
print(f"\nThe maximum number of chairs that can be occupied is {max_occupied}.")