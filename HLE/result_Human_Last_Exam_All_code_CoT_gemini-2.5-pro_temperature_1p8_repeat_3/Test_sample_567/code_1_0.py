def fibonacci(n):
    """Calculates the n-th Fibonacci number (F_1=1, F_2=1)."""
    if n <= 0:
        return 0
    elif n == 1:
        return 1
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b

# The values for 'a' where the embedding is volume-filling are of the form
# a = (F_{2k+1} / F_{2j+1})^2
# We are looking for the first non-trivial case. Let's take j=0 and k=1.
# The corresponding indices are 2*0+1=1 and 2*1+1=3.

j = 0
k = 1
idx1 = 2 * j + 1
idx2 = 2 * k + 1

F_idx1 = fibonacci(idx1)
F_idx2 = fibonacci(idx2)

# Calculate the value of a
a = (F_idx2 / F_idx1)**2

print(f"According to recent research, a value where the embedding obstruction becomes purely volume is given by the ratio of squared odd-indexed Fibonacci numbers.")
print(f"The first non-trivial value corresponds to j={j} and k={k}.")
print(f"The Fibonacci numbers are F_{idx1} = {F_idx1} and F_{idx2} = {F_idx2}.")
print(f"The equation for 'a' is a = (F_{idx2} / F_{idx1})^2.")
print(f"So, a = ({F_idx2} / {F_idx1})^2 = {a}")
print(f"Final calculated value for a: {int(a)}")
