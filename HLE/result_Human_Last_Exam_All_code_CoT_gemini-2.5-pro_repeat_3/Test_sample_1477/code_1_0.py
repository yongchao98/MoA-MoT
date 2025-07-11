import numpy as np
from scipy.optimize import minimize

def main():
    """
    This script demonstrates the calculation of the convergence bound for gradient
    descent on a logistic regression problem.
    """

    # 1. Define a sample non-separable dataset (x_i, y_i)
    # The data is not linearly separable if we can't find a 'w' such that
    # all y_i * w * x_i > 0.
    # For x=[1, -1] and y=[1, 1], we need w>0 and -w>0, which is a contradiction.
    x = np.array([1., 2., -1., -2.])
    y = np.array([1., -1., 1., -1.])
    N = len(x)
    print(f"Using sample data:\nx = {x}\ny = {y}\n")

    # 2. Define the required functions
    def sigma(t):
        """Sigmoid function."""
        # Clip t to avoid overflow in exp
        t = np.clip(t, -500, 500)
        return 1.0 / (1.0 + np.exp(-t))

    def R(w, x, y):
        """Loss function R(w)."""
        N = len(x)
        z = y * w * x
        # Add a small epsilon to log to avoid log(0)
        return -1/N * np.sum(np.log(sigma(z) + 1e-9))

    def R_prime(w, x, y):
        """First derivative of R(w)."""
        N = len(x)
        z = y * w * x
        return -1/N * np.sum(sigma(-z) * y * x)

    def R_double_prime(w, x, y):
        """Second derivative of R(w)."""
        N = len(x)
        z = y * w * x
        return 1/N * np.sum(sigma(z) * sigma(-z) * x**2)

    # 3. Calculate the uniform smoothness L
    # L is the maximum value of R''(w), which occurs at w=0.
    L = R_double_prime(0, x, y)
    
    # 4. Find the optimal w* by minimizing R(w)
    # We use a numerical optimizer for this.
    # The function to minimize needs to take a single vector argument.
    objective_func = lambda w_vec: R(w_vec[0], x, y)
    # The jacobian (gradient) for the optimizer.
    jacobian_func = lambda w_vec: np.array([R_prime(w_vec[0], x, y)])
    
    # Initial guess for w
    w_initial = np.array([0.0])
    
    result = minimize(objective_func, w_initial, jac=jacobian_func, method='BFGS')
    w_star = result.x[0]

    # 5. Calculate lambda, the smoothness at the optimum w*
    lambda_val = R_double_prime(w_star, x, y)
    
    # 6. Calculate the largest upper bound M for the learning rate gamma
    M = 2 / L

    # 7. Print the results
    print("--- Calculation Results ---")
    print(f"The uniform smoothness constant L is: {L:.4f}")
    print(f"The optimal weight w* is: {w_star:.4f}")
    print(f"The smoothness at the optimum, lambda, is: {lambda_val:.4f}")
    
    print("\n--- Convergence Condition ---")
    print("For gradient descent to be guaranteed to converge for any initialization,")
    print("the learning rate gamma must be less than M.")
    print("The equation for the largest upper bound M is:")
    print("M = 2 / L")
    print("\nFor the given data, the calculation is:")
    # Output each number in the final equation
    print(f"M = 2 / {L:.4f} = {M:.4f}")
    print("\nThis corresponds to answer choice C.")


if __name__ == "__main__":
    main()
