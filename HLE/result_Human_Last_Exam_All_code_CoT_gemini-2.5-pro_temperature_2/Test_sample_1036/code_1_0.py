# A brief explanation of the calculation:
# The problem reduces to counting the number of valid sets of remainders.
# A set of remainders is determined by which number 'x' from {0, ..., 99} is missing.
# We sum the number of valid assignment possibilities for each case of 'x'.

# Case x = 0: The set of remainders is {1, ..., 99}.
# This forces the assignment r_k = k-1, giving 1 solution.
count_x0 = 1

# Case x = 1: The number of valid assignments is 1.
count_x1 = 1

# Cases x = 2 to 98: The number of assignments for a given x is 2**(x-1).
# We sum these possibilities.
sum_x2_to_98 = 0
for x in range(2, 99):
    sum_x2_to_98 += 2**(x-1)

# Case x = 99: The number of valid assignments is 2**98.
count_x99 = 2**98

# The total number of solutions is the sum of counts from all cases.
# This sum simplifies to 2**99.
total_count = count_x0 + count_x1 + sum_x2_to_98 + count_x99

print(f"Number of solutions for case x=0: {count_x0}")
print(f"Number of solutions for case x=1: {count_x1}")
print(f"Number of solutions for cases x=2 to 98: {sum_x2_to_98}")
print(f"Number of solutions for case x=99: {count_x99}")
print("---")
print(f"{count_x0} + {count_x1} + {sum_x2_to_98} + {count_x99} = {total_count}")
print("---")
print(f"The total number of such positive integers is: {total_count}")