import numpy as np

def run_simulation():
    """
    Simulates SGD on a hard instance of stochastic logistic regression
    to empirically determine the convergence rate.
    """
    # We test the convergence for different numbers of iterations T.
    # We choose T such that the regime T = O(exp(D)) holds.
    # Let's set D = log(T), so T = exp(D).
    T_values = [1000, 4000, 16000, 64000]

    print("Running simulation to find the convergence rate.")
    print("We test if Error * T^alpha converges to a constant for alpha = 1/2, 2/3, or 1.\n")
    print(f"{'T':>8s} | {'D':>5s} | {'Error*T^(1/2)':>15s} | {'Error*T^(2/3)':>15s} | {'Error*T^1':>15s}")
    print("-" * 65)

    for T in T_values:
        # Set D according to the regime T = O(exp(D)). Let's use D=log(T).
        D = np.log(T)
        
        # Initial weight
        w = 0.0
        
        # Optimal weight for this problem is at the boundary
        w_star = -D
        
        # Loss function L(w) and its minimim value L(w*)
        loss_function = lambda v: np.log(1 + np.exp(v))
        min_loss = loss_function(w_star)

        # Step size schedule eta_t = c / sqrt(t)
        # The constant c can be tuned, but its value mainly affects the
        # constant factor of the rate, not the exponent. Let's use D
        # as a reasonable scaling factor for the learning rate, as suggested by theory.
        step_size_const = D 

        for t in range(1, T + 1):
            # Gradient of log(1 + exp(w)) is exp(w) / (1 + exp(w))
            gradient = np.exp(w) / (1 + np.exp(w))
            
            # Step size
            eta = step_size_const / np.sqrt(t)
            
            # SGD update
            w = w - eta * gradient
            
            # Projection onto the domain [-D, D]
            w = np.clip(w, -D, D)

        # Final excess loss (the error)
        final_loss = loss_function(w)
        error = final_loss - min_loss

        # Scale the error by powers of T to see which one becomes constant
        scaled_error_sqrt = error * (T**0.5)
        scaled_error_2_3 = error * (T**(2/3))
        scaled_error_1 = error * T

        print(f"{T:>8d} | {D:>5.2f} | {scaled_error_sqrt:>15.4f} | {scaled_error_2_3:>15.4f} | {scaled_error_1:>15.4f}")

run_simulation()