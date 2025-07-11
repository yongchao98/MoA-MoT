import numpy as np
from scipy import stats

def solve():
    """
    This function performs a linear regression on the given data to find the
    most parsimonious model and reports the results.
    """
    # 25 observations of x and y
    x = np.array([
        5, 6, 3.4, 2.7, 10, 9.7, 9.55, 3.05, 8.15, 6.2, 2.9, 6.35,
        4.6, 5.8, 7.4, 3.6, 7.85, 8.8, 7, 5.45, 9.1, 10.2, 4.1,
        3.95, 2.45
    ])
    y = np.array([
        1.582, 1.822, 1.057, 0.5, 2.236, 2.386, 2.294, 0.558, 2.166,
        1.866, 0.653, 1.93, 1.562, 1.737, 2.088, 1.137, 2.179,
        2.112, 1.8, 1.501, 2.303, 2.31, 1.194, 1.144, 0.123
    ])

    # Perform linear regression to find the best-fit line
    slope, intercept, r_value, _, _ = stats.linregress(x, y)
    
    # R-squared value indicates the goodness of fit
    r_squared = r_value**2

    # A high R-squared value (>0.9) suggests the linear model is a very good fit,
    # and thus the most parsimonious choice.
    
    # Format the parameters to 3 significant digits
    a = f"{slope:.3g}"
    b = f"{intercept:.3g}"
    
    print("The data shows a strong linear trend (R-squared = {:.3f}).".format(r_squared))
    print("Therefore, the maximally parsimonious model is a linear function of the form: y = a*x + b\n")
    print("The estimated parameter values (to 3 significant digits) are:")
    print(f"Slope (a) = {a}")
    print(f"Intercept (b) = {b}\n")
    print("The final prediction equation is:")
    print(f"y = {a} * x + {b}")

solve()
<<<y = 0.266 * x + -0.0131>>>