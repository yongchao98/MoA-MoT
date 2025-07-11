def fibonacci(n):
    """
    Computes the n-th Fibonacci number using the convention F_1=1, F_2=1.
    """
    if n <= 0:
        raise ValueError("Input must be a positive integer.")
    
    # Using the standard mathematical sequence F_1=1, F_2=1, F_3=2, ...
    # We can use a simple iterative approach.
    a, b = 1, 1 
    if n == 1:
        return a
    if n == 2:
        return b
    
    for _ in range(3, n + 1):
        a, b = b, a + b
    return b

# The value 'a' where the volume constraint ceases to be the only obstruction
# is the first corner of the McDuff-Schlenk staircase.
# This corner is located at a = (F_5 / F_3)^2, where F_n is the n-th Fibonacci number.

# Let's calculate the required Fibonacci numbers.
k = 1
numerator_index = 2 * k + 3  # which is 5
denominator_index = 2 * k + 1  # which is 3

F3 = fibonacci(denominator_index)
F5 = fibonacci(numerator_index)

# Now, calculate the value of 'a'.
a = (F5 / F3)**2

print("The problem is to find the value of 'a' where the symplectic embedding of the ellipsoid E(1,a)")
print("into a ball is no longer determined solely by the volume constraint.")
print("This transition occurs at the first corner of the McDuff-Schlenk 'infinite staircase'.")
print("The formula for this first corner is a = (F_5 / F_3)^2, where F_n are the Fibonacci numbers.")
print("")
print(f"The calculation is as follows:")
print(f"First, we find the required Fibonacci numbers (F_1=1, F_2=1, ...):")
print(f"F_3 = {F3}")
print(f"F_5 = {F5}")
print("")
print(f"Now, we compute the value of 'a':")
print(f"a = ({F5} / {F3})^2")
print(f"a = ({F5/F3})^2")
print(f"a = {a}")
