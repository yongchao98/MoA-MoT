import numpy as np
import matplotlib.pyplot as plt

def demonstrate_heavy_ball_convergence():
    """
    Analyzes the y-component dynamics of a proposed counterexample for the
    Heavy-ball method from Ghadimi et al. (2015). The analysis and simulation
    show that the y-component must converge to 0, contradicting the paper's
    claim of convergence to a non-stationary point with y=22.4.
    """
    # Parameters from Ghadimi et al. (2015) counterexample
    beta = 0.99
    gamma = 0.002

    # Initial conditions for the y-component.
    # We use non-zero initial conditions to demonstrate the convergence to 0.
    y_minus_1 = 10.0
    y_0 = 10.0

    # Number of iterations
    iterations = 3000

    # Store history for plotting
    y_history = [y_minus_1, y_0]

    # The update rule for the y-component is:
    # y_{k+1} = y_k + beta * (y_k - y_{k-1}) - gamma * y_k
    # y_{k+1} = (1 + beta - gamma) * y_k - beta * y_{k-1}
    
    y_k = y_0
    y_k_minus_1 = y_minus_1

    for i in range(iterations):
        y_k_plus_1 = (1 + beta - gamma) * y_k - beta * y_k_minus_1
        y_history.append(y_k_plus_1)
        
        # Update values for the next iteration
        y_k_minus_1 = y_k
        y_k = y_k_plus_1

    print("Demonstration of Heavy-ball method dynamics for a published counterexample.")
    print("The y-component of the gradient is simply 'y', so a stationary point must have y=0.")
    print(f"The paper claimed convergence to a point where y=22.4.")
    print(f"Our simulation shows the final value of y after {iterations} iterations is: {y_history[-1]:.6f}")
    print("This demonstrates that the y-component converges to 0, the stationary value.")
    
    # The plot will visually confirm the convergence to 0.
    # plt.figure(figsize=(10, 6))
    # plt.plot(y_history)
    # plt.title("Dynamics of the y-component in the Ghadimi et al. (2015) example")
    # plt.xlabel("Iteration k")
    # plt.ylabel("y_k")
    # plt.axhline(0, color='r', linestyle='--', label="Stationary value (y=0)")
    # plt.legend()
    # plt.grid(True)
    # plt.show()

demonstrate_heavy_ball_convergence()