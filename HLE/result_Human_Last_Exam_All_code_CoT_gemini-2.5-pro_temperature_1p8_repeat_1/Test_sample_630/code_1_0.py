import numpy as np

def sigmoid(z):
    """Numerically stable sigmoid function."""
    return np.where(z >= 0, 1 / (1 + np.exp(-z)), np.exp(z) / (1 + np.exp(z)))

def project_to_ball(w, D):
    """Projects a vector w onto the L2 ball of radius D."""
    norm = np.linalg.norm(w)
    if norm > D:
        return w * D / norm
    return w

def run_simulation():
    """
    Simulates SGD on a worst-case logistic regression problem and estimates
    the convergence rate.
    """
    # Problem parameters
    d = 2         # Dimension
    D = 5.0       # Radius of the parameter space W
    
    # Worst-case scenario: data x is always the same vector.
    # This makes the Hessian E[xx^T] singular.
    x_fixed = np.zeros(d)
    x_fixed[0] = 1.0

    # Loss function for this specific x: L(w) = log(1 + exp(w[0]))
    # The minimizer w_star is at the boundary, trying to make w[0] as small as possible.
    w_star = np.zeros(d)
    w_star[0] = -D
    loss_star = np.log(1 + np.exp(w_star[0]))

    # SGD parameters
    w_init = np.zeros(d) # Start at the origin
    
    # We will run SGD for various numbers of iterations (T)
    T_values = np.logspace(3, 6, 20).astype(int)
    excess_losses = []

    print("Running simulations for different T values...")
    for T in T_values:
        w = w_init.copy()
        # Use a learning rate schedule proportional to 1/sqrt(t)
        # eta_0 is a tuning parameter
        eta_0 = 1.0
        
        for t in range(1, T + 1):
            # Stochastic gradient (in our case, it's deterministic)
            grad = sigmoid(np.dot(w, x_fixed)) * x_fixed
            
            # Learning rate
            eta_t = eta_0 / np.sqrt(t)
            
            # SGD update step
            w = w - eta_t * grad
            
            # Projection step
            w = project_to_ball(w, D)
            
        # Calculate final excess loss
        final_loss = np.log(1 + np.exp(w[0]))
        excess_loss = final_loss - loss_star
        excess_losses.append(excess_loss)

    # Fit a line to the log-log plot of excess_loss vs. T
    # log(loss) = alpha * log(T) + const
    # We expect alpha to be close to -0.5, corresponding to a rate of 1/sqrt(T).
    log_T = np.log(T_values)
    log_loss = np.log(excess_losses)
    
    # Perform linear regression
    # A = np.vstack([log_T, np.ones(len(log_T))]).T
    # alpha, const = np.linalg.lstsq(A, log_loss, rcond=None)[0]
    # Using polyfit which is more common
    alpha, const = np.polyfit(log_T, log_loss, 1)

    print("\n--- Analysis of Convergence Rate ---")
    print(f"The relationship between excess loss and T is modeled as: Loss \u2248 C * T^\u03B1")
    print(f"This is equivalent to: log(Loss) \u2248 log(C) + \u03B1 * log(T)")
    print(f"We estimated the exponent \u03B1 by fitting a line to the log-log data.")
    print(f"Estimated slope (\u03B1): {alpha:.4f}")
    
    if abs(alpha + 0.5) < 0.1:
        print("\nThe estimated slope is close to -0.5, which supports a convergence rate of Theta(1/T^(1/2)).")
    elif abs(alpha + 1.0) < 0.1:
        print("\nThe estimated slope is close to -1.0, which supports a convergence rate of Theta(1/T).")
    else:
        print("\nThe estimated slope does not clearly match the common theoretical rates.")

run_simulation()
