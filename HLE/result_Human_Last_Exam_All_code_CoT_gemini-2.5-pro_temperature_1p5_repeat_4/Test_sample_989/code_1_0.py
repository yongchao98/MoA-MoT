import numpy as np

# This code demonstrates a counterexample to statement E.
# Statement E: Any strictly convex function has a unique global minimizer.
# The function f(x) = e^x is strictly convex, meaning its second derivative (e^x) is always positive.
# Let's test if this function has a global minimum by evaluating it for decreasing values of x.

print("This program will evaluate the strictly convex function f(x) = e^x at several points.")
print("We will show that as x gets smaller, f(x) also gets smaller, never reaching a minimum value.")
print("-" * 70)

# The following numbers are the inputs for our function f(x) = e^x
x_values = [10, 0, -10, -100, -1000]

# Here we evaluate the function and print the output
for x in x_values:
    # The equation is y = e^x
    y = np.exp(x)
    print(f"For the input x = {x}, the output of the equation y = e^x is: {y}")

print("-" * 70)
print("As you can see, as x approaches negative infinity, the function's value gets closer and closer to 0.")
print("However, for any finite x, e^x is always greater than 0.")
print("The function has a greatest lower bound (infimum) of 0, but it never reaches this value.")
print("Therefore, f(x) = e^x is a strictly convex function with no global minimizer.")
print("This demonstrates that statement E is false.")