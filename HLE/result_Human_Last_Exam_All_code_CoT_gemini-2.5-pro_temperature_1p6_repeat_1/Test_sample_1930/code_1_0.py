# The number of digits in the base. For the decimal system, this is 10.
b = 10

# The number of consecutive digits the function T(A)_n depends on.
# T(A)_n depends on A_n, A_{n+1}, A_{n+2}, so it's a window of 3 digits.
k = 3

# The formula for the dimension of the space of digitary functions is
# (b - 1) * b^(k - 1) + 1.
dimension = (b - 1) * b**(k - 1) + 1

# We print the equation with the values plugged in.
print(f"The dimension is calculated using the formula (b - 1) * b^(k - 1) + 1")
print(f"Dimension = ({b} - 1) * {b}^({k} - 1) + 1")
print(f"Dimension = {b-1} * {b}^{k-1} + 1")
print(f"Dimension = {b-1} * {b**(k-1)} + 1")
print(f"Dimension = {(b-1) * b**(k-1)} + 1")
print(f"Dimension = {dimension}")