def fibonacci(k):
    """Calculates the k-th Fibonacci number (F_1=1, F_2=1)."""
    if k <= 0:
        return 0
    a, b = 1, 1
    for _ in range(k - 2):
        a, b = b, a + b
    return b if k > 1 else 1

def number_of_ways(n):
    """
    Calculates the number of ways to tile the L-shaped region of size n.
    The formula is A_n = 2 * (F_{n-1})^2.
    """
    if n < 2:
        return 0
    fn_minus_1 = fibonacci(n - 1)
    return 2 * (fn_minus_1 ** 2)

# We present the answer as a formula using F_n notation.
# For example, to find the answer for n=4:
# n = 4
# F_n_minus_1_val = fibonacci(n - 1)
# result = number_of_ways(n)
# print(f"For n = {n}, the number of ways is {result}.")
# print(f"This is calculated as 2 * (F_{n-1})^2 = 2 * (F_3)^2 = 2 * {fibonacci(3)}^2 = {2*fibonacci(3)**2}")


print("The number of ways to tile the shape is given by the expression:")
print("2 * (F_{n-1})^2")
print("This can also be written as:")
print(f"{2} * F_{{n-{1}}} * F_{{n-{1}}}")
