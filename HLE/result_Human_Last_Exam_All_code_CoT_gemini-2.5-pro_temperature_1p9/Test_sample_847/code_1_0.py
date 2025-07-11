import math

def solve_cubic_root():
    """
    Solves the equation x^3 - x^2 - 1 = 0 for its real root using Newton's method.
    This equation arises from the cost analysis of the optimal strategy.
    """
    # The function f(x) = x^3 - x^2 - 1
    f = lambda x: x**3 - x**2 - 1
    # The derivative f'(x) = 3x^2 - 2x
    f_prime = lambda x: 3*x**2 - 2*x
    
    # Initial guess for the root. A quick check shows the root is between 1 and 2.
    x = 1.5
    
    # Iterate using Newton's method for high precision.
    for _ in range(20):
        x = x - f(x) / f_prime(x)
        
    return x

# Costs defined in the problem
cost_yes = 1
cost_no = 3
cost_compare = 2

# Strategy 1: Using only comparison questions.
# The asymptotic cost is (cost_compare / ln(2)) * n * ln(n).
# We calculate the coefficient of n*ln(n).
coeff_compare = cost_compare / math.log(2)

# Strategy 2: Using general yes/no questions.
# The analysis leads to the equation x^c_n = x^(c_n - c_y) + 1,
# which simplifies to x^3 = x^2 + 1 given the costs.
# First, find the real root of x^3 - x^2 - 1 = 0.
x_root = solve_cubic_root()

# The coefficient of the n*ln(n) term in the cost is c_y / ln(x).
coeff_general = cost_yes / math.log(x_root)

# The minimal cost is determined by the smaller of the two coefficients.
min_coeff = min(coeff_compare, coeff_general)

# To satisfy the instruction "output each number in the final equation",
# we print the cost parameters that define the problem and its solution.
print(f"Cost for a 'yes' answer: {cost_yes}")
print(f"Cost for a 'no' answer: {cost_no}")
print(f"Cost for a comparison: {cost_compare}")

# The problem asks for the minimal number of coins, which we interpret as the
# leading coefficient of the asymptotic cost function C(n) ~ A * n*ln(n).
# We output this coefficient A rounded to 3 decimal places.
print(f"The minimal cost coefficient is: {min_coeff:.3f}")
print("<<<2.616>>>")