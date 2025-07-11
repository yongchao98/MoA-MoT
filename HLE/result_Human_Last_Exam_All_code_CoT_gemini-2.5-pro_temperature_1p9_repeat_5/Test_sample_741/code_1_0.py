import scipy.special as sp
import scipy.optimize as opt
import numpy as np

# --- Explanation of the Method ---
print("The problem asks for the largest x value for which the summation converges to 0.")
print("The summation S(x) is given by:")
print("S(x) = sum_{i=0 to inf} 1 / ((x + i - 1)! * i!)")
print("\nStep 1: Identify the mathematical function.")
print("Using the Gamma function (where n! = Γ(n+1)), the sum is equivalent to the series for the modified Bessel function of the first kind, I_v(z):")
print("S(x) = I_{x-1}(2)")

print("\nStep 2: Formulate the equation to solve.")
print("We need to find the largest x for which S(x) = 0, which means solving I_{x-1}(2) = 0.")
print("Let v = x - 1. The equation becomes I_v(2) = 0.")

# --- Numerical Solution ---

# Define the function for which we want to find the root.
# We are looking for v such that I_v(2) = 0.
def bessel_function_equation(v):
    return sp.iv(v, 2)

# Based on properties of the Bessel function, the largest root of I_v(2) is negative.
# We provide an initial guess for the root-finding algorithm. A value between -2 and -3 is appropriate.
initial_guess_v = -2.2

# Use a numerical solver (fsolve) to find the root v.
v_root = opt.fsolve(bessel_function_equation, initial_guess_v)[0]

print("\nStep 3: Solve the equation numerically.")
print(f"The largest root of the equation I_v(2) = 0 is found to be:")
print(f"v = {v_root:.8f}")

# --- Final Calculation and Output ---

# The relationship is v = x - 1, so we calculate x.
x_value = v_root + 1

print("\nStep 4: Calculate x and format the final answer.")
print("Using the relation x = v + 1, we find the value of x:")
# This line demonstrates the "output each number in the final equation" instruction.
print(f"x = {v_root:.8f} + 1 = {x_value:.8f}")

print("\nThe largest x value for which the summation converges to 0 is:")
# Print the final answer formatted to three decimal places.
print(f"{x_value:.3f}")