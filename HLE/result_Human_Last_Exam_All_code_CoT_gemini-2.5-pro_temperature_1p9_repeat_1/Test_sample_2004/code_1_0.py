import numpy as np
from scipy import stats

def solve():
    """
    This function finds the maximally parsimonious model for the given data,
    calculates the parameters, and prints the final equation.
    """
    # The 25 observations of x and y
    x_data = np.array([5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35, 4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1, 3.95, 2.45])
    y_data = np.array([1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166, 1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179, 2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123])

    # Based on analysis, the best-fit model has the form: y = a * sqrt(x) + b
    # We linearize this by using z = sqrt(x), so y = a * z + b
    z_data = np.sqrt(x_data)

    # Perform linear regression to find 'a' (slope) and 'b' (intercept)
    # This minimizes the sum of squared differences.
    slope, intercept, r_value, p_value, std_err = stats.linregress(z_data, y_data)
    
    # Format the parameters to 3 significant digits
    # 'g' format specifier is used for general format, which works well for significant figures
    a_est = f"{slope:.3g}"
    b_val = intercept
    
    # Build the final equation string, handling the sign of the intercept
    if b_val >= 0:
        sign = "+"
        b_est = f"{b_val:.3g}"
    else:
        sign = "-"
        # Use abs() for the number part of the string
        b_est = f"{abs(b_val):.3g}"

    print(f"The maximally parsimonious model for the data is of the form y = a*sqrt(x) + b.")
    print("After performing a least-squares fit, the estimated equation is:")
    print(f"y = {a_est} * sqrt(x) {sign} {b_est}")

solve()
<<<y = 0.963 * sqrt(x) - 0.585>>>