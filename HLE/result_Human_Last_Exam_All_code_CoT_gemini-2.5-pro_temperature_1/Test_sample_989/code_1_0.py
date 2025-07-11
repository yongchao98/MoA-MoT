import numpy as np
from scipy.optimize import minimize

# Statement E is: "Any strictly convex function has a unique global minimizer"
# We will show this is false with a counterexample.
# The function f(x) = e^x is a strictly convex function because its second 
# derivative, f''(x) = e^x, is positive for all real x.
# However, this function does not have a global minimum. It approaches 0 as x -> -infinity, but never reaches it.
# We can demonstrate this by asking a numerical optimizer to find the minimum.

print("Testing Statement E: 'Any strictly convex function has a unique global minimizer'.")
print("We will use the strictly convex function f(x) = exp(x) as a counterexample.")
print("This function's infimum (greatest lower bound) is 0, but it never reaches a minimum value.")
print("Let's see what happens when we ask a numerical optimizer to find the minimum, starting from x=10.")
print("-" * 50)

# Define the strictly convex function
def f(x):
    return np.exp(x[0])

# An initial guess for the minimizer
x0 = [10.0]

# Attempt to find the minimum using scipy.optimize.minimize
result = minimize(f, x0, method='BFGS')

print("Optimizer Result Details:")
print(result)
print("-" * 50)

print("Analysis of the Result:")
print(f"The optimization terminated. Was it successful? {result.success}")
print(f"The reason for termination is: '{result.message}'")
print(f"Final value of x: {result.x[0]:.2f}")
print(f"Final function value f(x): {result.fun:.2e}")
print("\nConclusion:")
print("The optimizer did not find a true minimum. It stopped because the gradient became very small, a consequence of the function value getting extremely close to its limit of 0.")
print(f"The final x value is a large negative number ({result.x[0]:.2f}), and the function value ({result.fun:.2e}) is nearly zero.")
print("This shows that for f(x) = e^x, there is no value x* for which f(x*) is the minimum. The function decreases indefinitely.")
print("Therefore, the statement that *any* strictly convex function *has* a global minimizer is false.")
