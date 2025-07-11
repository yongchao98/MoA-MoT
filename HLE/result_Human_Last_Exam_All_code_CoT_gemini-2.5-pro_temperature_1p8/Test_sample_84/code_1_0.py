import math

def calculate_alpha(power1, power2):
    """
    Calculates the exponent alpha in the asymptotic growth rate of the polynomial degree.
    The degree d_n is proportional to the square root of the product of the lengths of the intervals.
    The first interval has length proportional to n^power1.
    The second interval has length proportional to n^power2.
    d_n = Theta(n^alpha)
    """
    
    # The lengths of the intervals are approximately n^2 and n^10.
    # L1 ~ n^power1
    # L2 ~ n^power2
    # d_n ~ sqrt(L1 * L2) = sqrt(n^power1 * n^power2) = sqrt(n^(power1 + power2))
    # d_n ~ n^((power1 + power2) / 2)
    
    alpha = (power1 + power2) / 2.0
    return alpha

# The sets of integers are {1, 2, ..., n^2} and {n^2+1, ..., n^10}.
# The length of the first interval is of order n^2.
power1 = 2
# The length of the second interval is of order n^10.
power2 = 10

# Calculate alpha
alpha = calculate_alpha(power1, power2)

# Print the explanation and the final equation
print(f"The length of the first set of points is of the order n^{power1}.")
print(f"The length of the second set of points is of the order n^{power2}.")
print("The degree d_n of the polynomial behaves as the square root of the product of the lengths of these intervals.")
print(f"d_n = Theta(sqrt(n^{power1} * n^{power2})) = Theta(sqrt(n^({power1}+{power2}))) = Theta(n^({(power1+power2)}/2))")
print(f"Thus, alpha = ({power1} + {power2}) / 2 = {alpha}")
print("The determined value of alpha is:")
print(alpha)
