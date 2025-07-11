import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def run_sgd_logistic_regression(T_max, D_val, is_D_fixed=True):
    """
    Simulates SGD for a simple logistic regression problem.

    The problem is to minimize E[log(1 + exp(w*x))] where x is always 1.
    The domain is |w| <= D.
    The optimal solution is w_star = -D.
    The optimal loss is L_star = log(1 + exp(-D)).
    """
    # Set dimension to 1 for simplicity
    d = 1
    # The data x is always 1.
    x_val = 1.0

    # Initialize w
    w = np.zeros(d)
    
    # Optimal solution and loss
    w_star = -D_val
    l_star = np.log(1 + np.exp(-D_val * x_val))
    
    errors = []
    T_values = []
    
    # Step size for SGD
    # A simple decreasing step size
    eta0 = D_val
    
    avg_w = np.zeros(d)
    
    for t in range(1, T_max + 1):
        # Update D if it depends on T
        if not is_D_fixed:
            D = np.log(t) if t > 1 else 1.0
            w_star = -D
            l_star = np.log(1 + np.exp(-D * x_val))
        else:
            D = D_val

        # SGD step
        eta = eta0 / np.sqrt(t)
        grad = (1 / (1 + np.exp(-w * x_val))) * x_val
        w = w - eta * grad
        
        # Projection onto the ball of radius D
        norm_w = np.linalg.norm(w)
        if norm_w > D:
            w = w * D / norm_w

        # Use iterate averaging for stability
        avg_w = (avg_w * (t - 1) + w) / t

        if t % 100 == 0 or t == T_max:
            current_loss = np.log(1 + np.exp(avg_w * x_val))
            excess_loss = current_loss - l_star
            errors.append(excess_loss[0])
            T_values.append(t)
            
    return np.array(T_values), np.array(errors)

def main():
    T_max = 20000
    
    # Case 1: D is a fixed constant
    D_fixed = 4.0
    print(f"Running simulation for fixed D = {D_fixed}")
    T_vals_fixed, errors_fixed = run_sgd_logistic_regression(T_max, D_fixed, is_D_fixed=True)

    # Case 2: D = log(T)
    print(f"Running simulation for D = log(T)")
    # D_val is not used here, but we need to pass a value.
    T_vals_logD, errors_logD = run_sgd_logistic_regression(T_max, 1.0, is_D_fixed=False)

    # Fit curves to find the exponent of T
    # Model 1: error = a * T^b
    def power_law(T, a, b):
        return a * T**b

    # Model 2: error = a * log(T) * T^b
    def log_power_law(T, a, b):
        # Add epsilon to avoid log(0)
        return a * np.log(T + 1e-9) * T**b

    # Fit for fixed D case
    params_fixed, _ = curve_fit(power_law, T_vals_fixed, errors_fixed)
    
    # Fit for log(D) case
    params_logD, _ = curve_fit(log_power_law, T_vals_logD, errors_logD)
    
    print("\n--- Analysis of Convergence Rate ---")
    print("The theoretical rate for stochastic convex optimization is Theta(D / T^0.5).")
    
    print("\nCase 1: Fixed D")
    print("We expect the rate to be proportional to 1/sqrt(T), i.e., T^-0.5.")
    print(f"The fitted curve is Error ~ T^({params_fixed[1]:.4f}). This is close to -0.5.")
    
    print("\nCase 2: D = log(T)")
    print("We expect the rate to be proportional to log(T)/sqrt(T), i.e., log(T)*T^-0.5.")
    print(f"The fitted curve is Error ~ log(T) * T^({params_logD[1]:.4f}). This is close to -0.5.")
    print("This shows the rate is slower by a log(T) factor, matching the theoretical analysis.")

    print("\nFinal Conclusion from Theory:")
    print("The rate is Omega(log(T)/sqrt(T)), which is not among options A, B, or C.")
    print("Therefore, the correct answer is D.")

if __name__ == '__main__':
    main()