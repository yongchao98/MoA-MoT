def double_factorial(n):
    """
    Calculates the double factorial (n!!).
    Example: 5!! = 5 * 3 * 1 = 15
    """
    if n < 0:
        return 0
    if n == 0 or n == -1:
        return 1
    result = 1
    for i in range(n, 0, -2):
        result *= i
    return result

# The task is to find a(3). We need a(1) and a(2) first.
# The recurrence relation is: a(n) = (2n-1)!! - sum_{k=1}^{n-1} a(k) * (2(n-k)-1)!!

# 1. Calculate a(1)
# a(1) = (2*1-1)!! = 1!!
a1 = double_factorial(1)

# 2. Calculate a(2)
# a(2) = (2*2-1)!! - a(1)*(2*(2-1)-1)!! = 3!! - a(1)*1!!
df_3 = double_factorial(3)
df_1 = double_factorial(1)
a2 = df_3 - (a1 * df_1)

# 3. Calculate a(3)
# a(3) = (2*3-1)!! - (a(1)*(2*(3-1)-1)!! + a(2)*(2*(3-2)-1)!!)
# a(3) = 5!! - (a(1)*3!! + a(2)*1!!)
df_5 = double_factorial(5)
sum_term = (a1 * df_3) + (a2 * df_1)
a3 = df_5 - sum_term

# Print the final equation with all numbers substituted, as requested.
print("To find a(3), we calculate the following equation:")
print(f"a(3) = {df_5} - ({a1} * {df_3} + {a2} * {df_1})")
print(f"The result is: {a3}")