import math

def f(x):
    """A simple linear function."""
    return x

def grad_f(x):
    """The gradient of f(x)."""
    return 1

def gradient_descent_demo():
    """
    Demonstrates that Gradient Descent can converge to a non-stationary point
    if the step-size sequence is summable.
    """
    # Initial point
    x0 = 5.0
    xk = x0
    
    # Number of iterations
    iterations = 100000
    
    print(f"Running Gradient Descent for f(x)=x, starting at x0 = {x0}")
    print("Using a summable step-size sequence gamma_k = 1/k^2 for k>=1.")
    print("-" * 50)

    # Perform iterations
    for k in range(1, iterations + 1):
        gamma_k = 1.0 / (k**2)
        gradient = grad_f(xk)
        xk = xk - gamma_k * gradient

    # The theoretical limit point
    # The sum of 1/k^2 from k=1 to infinity is pi^2 / 6
    x_star_theoretical = x0 - (math.pi**2) / 6
    
    # The calculated limit point
    x_star_actual = xk
    
    # The gradient at the limit point
    grad_at_limit = grad_f(x_star_actual)

    print(f"Algorithm has run for {iterations} iterations.")
    print(f"Final point x_k = {x_star_actual}")
    print(f"Theoretical limit point x* = {x_star_theoretical}")
    print(f"The sequence converges, and the difference is: {abs(x_star_actual - x_star_theoretical)}")
    print("-" * 50)
    print("Now, let's check for stationarity at the limit point.")
    print(f"The gradient at the limit point is nabla f(x*) = {grad_at_limit}")
    
    if grad_at_limit != 0:
        print("Since the gradient is non-zero, the limit point is NOT a first-order stationary point.")
    else:
        print("The gradient is zero, so the limit point is a stationary point.")

gradient_descent_demo()
