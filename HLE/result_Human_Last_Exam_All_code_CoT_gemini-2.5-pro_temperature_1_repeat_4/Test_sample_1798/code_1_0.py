import numpy as np
from scipy.integrate import quad

# Step 1 & 2: Define the function f(x) based on the circle equation
# The circle passes through A(1,0), B(10,0), C(4,3).
# Center (h,k), radius r.
# Perpendicular bisector of AB is x = (1+10)/2 = 5.5, so h=5.5.
# (1-5.5)^2 + (0-k)^2 = r^2  => 20.25 + k^2 = r^2
# (4-5.5)^2 + (3-k)^2 = r^2  => 2.25 + 9 - 6k + k^2 = r^2
# 20.25 = 11.25 - 6k => 9 = -6k => k = -1.5
# r^2 = 20.25 + (-1.5)^2 = 20.25 + 2.25 = 22.5
h = 5.5
k = -1.5
r_squared = 22.5

# The equation of the circle is (x - 5.5)^2 + (y + 1.5)^2 = 22.5
# Solving for y gives y = -1.5 +/- sqrt(22.5 - (x-5.5)^2).
# Since it must pass through C(4,3), we take the positive root.
# f(x) = -1.5 + sqrt(22.5 - (x-5.5)^2)
def f(x):
  return -1.5 + np.sqrt(22.5 - (x - 5.5)**2)

# Step 3: Find the normalization constant alpha
# alpha = 1 / integral of f(x) from 1 to 10.
total_integral, _ = quad(f, 1, 10)
alpha = 1 / total_integral

# Step 4: Calculate the probability P(X < 3)
# P(X<3) = alpha * integral of f(x) from 1 to 3
prob_integral, _ = quad(f, 1, 3)
probability = alpha * prob_integral

# Output the results
print("Step 1: Find the equation of the circle.")
print(f"The center of the circle is (h, k) = ({h}, {k}).")
print(f"The radius squared is r^2 = {r_squared}.")
print(f"The equation of the circle is (x - {h})^2 + (y - {-k})^2 = {r_squared}\n")

print("Step 2: Find the normalization constant alpha.")
print(f"The normalization condition is: alpha * integral(f(x) from 1 to 10) = 1")
print(f"The value of the integral from 1 to 10 is: {total_integral}")
print(f"alpha = 1 / {total_integral} = {alpha}\n")

print("Step 3: Calculate the probability P(X < 3).")
print(f"P(X < 3) = alpha * integral(f(x) from 1 to 3)")
print(f"The value of the integral from 1 to 3 is: {prob_integral}")
print(f"P(X < 3) = {alpha} * {prob_integral} = {probability}\n")

print("Final Answer:")
print(f"The value of alpha is: {alpha}")
print(f"The value of P(X < 3) is: {probability}")
<<<alpha: 0.04683070113063544, P(X < 3): 0.14716738906935742>>>