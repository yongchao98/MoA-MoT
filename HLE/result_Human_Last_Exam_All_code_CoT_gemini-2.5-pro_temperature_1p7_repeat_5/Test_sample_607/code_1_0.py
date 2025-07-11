import numpy as np
import matplotlib.pyplot as plt

def f(x):
    """
    A smooth function designed to have a "trap" for the heavy-ball method.
    It's a quadratic function plus a sine term.
    f(x) = 0.05 * x^2 - cos(10*x)
    Stationary points are solutions to f'(x) = 0.1*x + 10*sin(10*x) = 0.
    One such point is x=0.
    The point x = pi/10 (approx 0.314) is not a stationary point.
    """
    return 0.05 * x**2 - np.cos(10 * x)

def df(x):
    """
    Derivative of f(x).
    """
    return 0.1 * x + 10 * np.sin(10 * x)

def heavy_ball_method(x_start, x_prev, gamma, beta, num_steps):
    """
    Implements the Heavy-ball method.
    """
    x_current = x_start
    x_previous = x_prev
    
    # Store history for plotting
    history = [x_previous, x_current]
    
    print(f"Heavy-ball method started at x_0={x_start}, x_(-1)={x_prev}")
    print("Targeting convergence to a non-stationary point.\n")
    
    for k in range(num_steps):
        # Heavy-ball update rule
        gradient = df(x_current)
        x_next = x_current + beta * (x_current - x_previous) - gamma * gradient
        
        # Update for next iteration
        x_previous = x_current
        x_current = x_next
        history.append(x_current)

    final_x = history[-1]
    final_grad = df(final_x)
    
    print(f"Final point after {num_steps} steps: x* = {final_x}")
    print(f"Gradient at the final point: ∇f(x*) = {final_grad}")
    
    if np.abs(final_grad) > 1e-3:
        print("The algorithm converged to a point that is NOT first-order stationary.")
    else:
        print("The algorithm converged to a stationary point.")

    # Plotting the results
    x_vals = np.linspace(-0.5, 0.5, 400)
    y_vals = f(x_vals)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, label="f(x)")
    plt.plot(history, f(np.array(history)), 'r-o', label="HB Path", markersize=4)
    plt.title("Heavy-Ball Method Convergence")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.axvline(np.pi / 10, color='g', linestyle='--', label=f'Non-stationary point x={np.pi/10:.3f}')
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    # A point that is not stationary. For example, x=pi/10 (~0.314).
    # Here, f'(pi/10) = 0.1*pi/10 + 10*sin(pi) = pi/100 != 0
    non_stationary_point = np.pi / 10

    # These parameters and starting points are carefully chosen
    # to make the Heavy-ball method converge near the non-stationary point.
    gamma = 0.095
    beta = 0.9
    num_steps = 100
    
    # Start close to the non-stationary point
    x_start = non_stationary_point + 0.001
    x_prev = non_stationary_point - 0.001
    
    heavy_ball_method(x_start, x_prev, gamma, beta, num_steps)