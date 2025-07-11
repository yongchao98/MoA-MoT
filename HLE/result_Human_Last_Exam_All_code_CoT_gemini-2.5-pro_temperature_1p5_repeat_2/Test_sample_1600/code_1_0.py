import collections

# Use a dictionary for memoization to store the results of a(n)
memo = {0: 1, 1: 1}

def a(n):
    """
    Calculates the number of non-vanishing Feynman diagrams of order 2n
    for electron or photon propagators in QED using a recurrence relation.
    a(n) = (2*n - 3) * a(n-1) + (n - 1)**2 * a(n-2) for n >= 2
    Base cases: a(0) = 1, a(1) = 1.
    """
    if n in memo:
        return memo[n]
    
    # Recursive step
    val = (2 * n - 3) * a(n - 1) + (n - 1)**2 * a(n - 2)
    memo[n] = val
    return val

# The user wants to find a(3)
n_val = 3
result = a(n_val)

# To show the final equation, we need the values of a(1) and a(2)
a_2 = memo[2]
a_1 = memo[1]

# Output the final equation with the computed values
print(f"The number of diagrams a(n) can be found using the recurrence relation:")
print(f"a(n) = (2*n - 3) * a(n-1) + (n - 1)^2 * a(n-2)")
print(f"\nFor n = {n_val}, the equation is:")
print(f"a({n_val}) = (2*{n_val} - 3) * a({n_val-1}) + ({n_val} - 1)^2 * a({n_val-2})")
print(f"a({n_val}) = ({2*n_val - 3}) * a({n_val-1}) + ({n_val - 1}^2) * a({n_val-2})")
print(f"a({n_val}) = ({2*n_val - 3}) * {a_2} + ({n_val - 1}**2) * {a_1}")
print(f"a({n_val}) = { (2*n_val - 3) * a_2 } + { (n_val - 1)**2 * a_1 }")
print(f"a({n_val}) = {result}")
