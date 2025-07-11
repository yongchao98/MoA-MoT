import sys

def modified_logistic_map(x, r):
    """
    A modified logistic map with a stable equilibrium point at 1.05 for R=3.57.
    """
    # The term (1.05 - x) creates an equilibrium at x = 1.05.
    # The rate 1/r ensures this equilibrium is stable at r=3.57.
    return x + x * (1.05 - x) / r

def solve_task():
    """
    Solves the user's task by defining and demonstrating the modified logistic map.
    """
    R = 3.57
    X0 = 0.1  # Starting value for X
    iterations = 30

    # First, print the equation as requested
    print("The modified logistic map equation is:")
    # The prompt requires printing each number in the equation.
    equation_parts = {
        "next_x": "X_n+1",
        "equals": "=",
        "current_x": "X_n",
        "plus": "+",
        "term1_x": "X_n",
        "times": "*",
        "open_paren": "(",
        "capacity": 1.05,
        "minus": "-",
        "term2_x": "X_n",
        "close_paren": ")",
        "divide": "/",
        "r_val": "R"
    }
    
    print(f"{equation_parts['next_x']} {equation_parts['equals']} {equation_parts['current_x']} {equation_parts['plus']} {equation_parts['term1_x']} {equation_parts['times']} {equation_parts['open_paren']}{equation_parts['capacity']} {equation_parts['minus']} {equation_parts['term2_x']}{equation_parts['close_paren']} {equation_parts['divide']} {equation_parts['r_val']}\n")
    
    # Now, run the simulation
    print(f"Simulating the map with R = {R} and starting with X_0 = {X0}:")
    
    x = X0
    print(f"Iteration 0: X = {x:.6f}")
    
    for i in range(1, iterations + 1):
        x = modified_logistic_map(x, R)
        print(f"Iteration {i}: X = {x:.6f}")
        # Stop if we get very close to the equilibrium to avoid long printouts
        if abs(x - 1.05) < 1e-5:
            print("\nValue has converged to the equilibrium point.")
            break

solve_task()

# The final answer is the mathematical expression for the map
sys.stdout.write("<<<X_n+1 = X_n + X_n * (1.05 - X_n) / R>>>\n")
