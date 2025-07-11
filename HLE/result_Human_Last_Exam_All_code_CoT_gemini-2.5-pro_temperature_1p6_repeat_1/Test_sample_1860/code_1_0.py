import math

# Based on the analytical solution, the equation f(x) = g(x) has 8 roots
# when k is in the interval [1/3, sqrt(2)/4).
# This script calculates and prints the bounds of this interval.

# Define the components of the lower bound
k_lower_num = 1
k_lower_den = 3
k_lower_val = k_lower_num / k_lower_den

# Define the components of the upper bound
k_upper_sqrt_arg = 2
k_upper_den = 4
k_upper_val = math.sqrt(k_upper_sqrt_arg) / k_upper_den

print("The range of values for k is determined by analyzing the number of intersections between the graphs of f(x) and g(x).")
print("The analysis identifies two critical values for k that change the number of roots.")

print("\nThe lower bound of the range for k is inclusive.")
print(f"The equation for the lower bound is k = {k_lower_num} / {k_lower_den}")
print(f"This corresponds to the case where the line g(x) intersects the endpoints of the semi-circles f(x).")

print("\nThe upper bound of the range for k is exclusive.")
print(f"The equation for the upper bound is k = sqrt({k_upper_sqrt_arg}) / {k_upper_den}")
print(f"This corresponds to the case where the line g(x) is tangent to the semi-circles f(x).")

print("\nFinal Result:")
print(f"The range of values for k is [{k_lower_num}/{k_lower_den}, sqrt({k_upper_sqrt_arg})/{k_upper_den}).")