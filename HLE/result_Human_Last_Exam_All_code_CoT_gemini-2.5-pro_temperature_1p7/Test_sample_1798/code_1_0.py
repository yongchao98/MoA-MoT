import math
from scipy.integrate import quad

# Step 1 & 2: Define the function f(x) from the circle's equation
# The circle passes through A(1,0), B(10,0), and C(4,3).
# Center's x-coordinate h = (1+10)/2 = 5.5
# Plugging A(1,0) into (x-h)^2 + (y-k)^2 = r^2 gives:
# (1-5.5)^2 + (0-k)^2 = r^2 => 20.25 + k^2 = r^2
# Plugging C(4,3) into (x-h)^2 + (y-k)^2 = r^2 gives:
# (4-5.5)^2 + (3-k)^2 = r^2 => 2.25 + 9 - 6k + k^2 = r^2
# Equating the two expressions for r^2:
# 20.25 + k^2 = 11.25 - 6k + k^2 => 9 = -6k => k = -1.5
# Solving for r^2: r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5
# The circle equation is (x - 5.5)^2 + (y + 1.5)^2 = 22.5
# To get y = f(x), we solve for y: y = -1.5 +/- sqrt(22.5 - (x-5.5)^2)
# Checking point C(4,3): 3 = -1.5 + sqrt(22.5 - (4-5.5)^2) = -1.5 + sqrt(20.25) = -1.5 + 4.5. This is correct.
# So we use the positive branch.
h = 5.5
k = -1.5
r_squared = 22.5

def f(x):
    """The equation for the upper arc of the circle."""
    return k + math.sqrt(r_squared - (x - h)**2)

# Step 3: Calculate the normalization constant alpha
# We need to solve Integral[alpha * f(x), {x, 1, 10}] = 1
# This means alpha = 1 / Integral[f(x), {x, 1, 10}]
try:
    total_integral, _ = quad(f, 1, 10)
    alpha = 1 / total_integral

    # Step 4: Calculate the probability P(X<3)
    # P(X<3) = Integral[alpha * f(x), {x, 1, 3}] = alpha * Integral[f(x), {x, 1, 3}]
    # The domain of X is [1, 10], so P(X<3) corresponds to the integral from 1 to 3.
    partial_integral, _ = quad(f, 1, 3)
    prob_less_than_3 = alpha * partial_integral

    # Print the results
    print("The equation of the circle is (x - h)^2 + (y - k)^2 = r^2")
    print(f"h = {h}")
    print(f"k = {k}")
    print(f"r^2 = {r_squared}")
    print("\nCalculated values:")
    print(f"alpha = {alpha}")
    print(f"P(X<3) = {prob_less_than_3}")

except ImportError:
    print("This script requires the scipy library.")
    print("Please install it using: pip install scipy")
except Exception as e:
    print(f"An error occurred: {e}")
