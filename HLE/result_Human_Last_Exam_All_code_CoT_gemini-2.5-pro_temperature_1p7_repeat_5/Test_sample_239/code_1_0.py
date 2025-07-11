import math

# Question 1: How many lines of code have compiler errors?
# 1. `unsigned long long...` -> long long is not a valid XVM type.
# 2. `scanf("%d %d %d", ...)` -> %d is the wrong format specifier for large numbers.
# 3. `printf("%d", ...)` -> %d is the wrong format specifier for the large result.
q1_answer = 3

# Question 2: What is the total memory (in D) used for variables in an optimized rewrite?
# The minimal set of variables is n, m, and a to store the input.
# To hold values up to 4*10^9, they must be of type `unsigned long`.
# Size of `unsigned long` = 12D.
# Total memory = 3 variables * 12D/variable.
q2_answer = 3 * 12

# Question 3: Can this problem be solved without branching statements?
# Yes, by using ceiling division formula: (n + a - 1) // a in integer arithmetic.
q3_answer = "Y"

# Question 4: What is the expected output for m = n = 4*10^9 and a = 1?
n = 4 * 10**9
m = 4 * 10**9
a = 1

# The number of flagstones for a side is the ceiling of the length divided by the flagstone size.
# In Python, we can use math.ceil for clarity, or integer division trick `(length + size - 1) // size`.
na = math.ceil(n / a)
ma = math.ceil(m / a)

# The result is the product of the number of flagstones for each dimension.
# Python's integers can handle arbitrary size, so overflow is not an issue here.
result = na * ma

# Format the final part of the answer string as an equation.
q4_answer = f"{int(na)} * {int(ma)} = {int(result)}"

# Print the final answer in the format x:y:z:t
print(f"{q1_answer}:{q2_answer}:{q3_answer}:{q4_answer}")
