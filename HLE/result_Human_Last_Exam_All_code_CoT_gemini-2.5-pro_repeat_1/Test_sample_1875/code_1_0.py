import numpy as np

def solve_regression():
    """
    Constructs a parsimonious model to predict y from x1, x2, and x3.
    """
    # 1. Define the data from the problem description.
    data = np.array([
        [-1, -1, -1, 34.3],
        [ 1, -1, -1, 94.6],
        [-1,  1, -1, 52.5],
        [ 1,  1, -1, 75.1],
        [-1, -1,  1, 28.4],
        [ 1, -1,  1, 48.3],
        [-1,  1,  1,  7.0],
        [ 1,  1,  1, 70.8]
    ])

    y = data[:, 3]
    x1 = data[:, 0]
    x2 = data[:, 1]
    x3 = data[:, 2]
    
    N = len(y)

    # 2. Construct the full model design matrix X.
    # The columns are: Intercept, x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3
    X = np.c_[np.ones(N), x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3]
    term_names = ["Intercept", "x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3"]

    # 3. Calculate the coefficients for the full model.
    # Due to orthogonality, beta = (X'X)^-1 * X'y = (1/N) * X'y
    coefficients = (1/N) * X.T @ y

    # 4. Assess effects and select a parsimonious model.
    # Let's look at the calculated coefficients:
    # b0 (Intercept): 51.375
    # b1 (x1):        20.825
    # b2 (x2):        -0.025
    # b3 (x3):       -12.625
    # b12 (x1*x2):     0.775
    # b13 (x1*x3):   -17.600
    # b23 (x2*x3):     2.175
    # b123 (x1*x2*x3): 6.450
    #
    # The coefficients for x2, x1*x2, and x2*x3 are very small compared to others.
    # The largest effects are for x1, x3, and their interaction x1*x3.
    # The three-way interaction x1*x2*x3 is moderately large, but including it
    # without its parent main effect x2 would violate the hierarchy principle.
    # A strong, simple, and hierarchical model includes the intercept and the three largest effects.
    
    b0 = coefficients[0]
    b1 = coefficients[1]
    b3 = coefficients[3]
    b13 = coefficients[5]

    # 5. Format the final equation string.
    # The model is y = b0 + b1*x1 + b3*x3 + b13*x1*x3
    
    # We will build the string piece by piece to handle signs correctly.
    # Round coefficients to one decimal place.
    b0_r = round(b0, 1)
    b1_r = round(b1, 1)
    b3_r = round(b3, 1)
    b13_r = round(b13, 1)

    equation = f"y = {b0_r}"
    equation += f" + {b1_r}*x1" if b1_r >= 0 else f" - {abs(b1_r)}*x1"
    equation += f" + {b3_r}*x3" if b3_r >= 0 else f" - {abs(b3_r)}*x3"
    equation += f" + {b13_r}*x1*x3" if b13_r >= 0 else f" - {abs(b13_r)}*x1*x3"
    
    print(equation)

solve_regression()