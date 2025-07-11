import math

def solve():
    """
    This function modifies the logistic map to create a stable equilibrium point
    at approximately 1.05 for R=3.57, and demonstrates its convergence.
    """
    R = 3.57
    # Initial value for X, can be almost any value between 0 and 2
    x_current = 0.5

    # The modified logistic map function.
    # The term a = R**(2/3) is chosen to precisely place the equilibrium
    # point at approximately 1.05 for the given R=3.57.
    def modified_logistic_map(x, r):
        a = r**(2.0/3.0)
        return r * x / (1 + a * x**2)

    print("--- Modified Logistic Map Simulation ---")
    print("\nOriginal Map (for context): X_n+1 = R * X_n * (1 - X_n)")
    print("At R=3.57, the original map is chaotic.")

    print("\nModified Map Equation:")
    print("X_n+1 = R * X_n / (1 + R**(2/3) * X_n**2)")

    print(f"\nSimulating with R = {R}...")
    print("-" * 40)
    print(f"Iteration 0:  X = {x_current:.6f}")

    # Iterate to show the system converging to the stable point
    for n in range(1, 21):
        x_current = modified_logistic_map(x_current, R)
        print(f"Iteration {n:<2}: X = {x_current:.6f}")

    print("-" * 40)

    # Calculate the theoretical equilibrium point from the formula X^2 = (R-1)/R^(2/3)
    eq_val = math.sqrt((R - 1) / (R**(2.0/3.0)))
    
    print(f"The system converges to a stable equilibrium point at X ≈ {eq_val:.6f}.")
    print("This value is approximately equal to the target of 1.05.")

    print("\n--- Final Equation with Numbers ---")
    print("The equilibrium is reached when X_n+1 = X_n.")
    print("Let X_eq be the equilibrium value.")
    print("X_eq = R * X_eq / (1 + R**(2/3) * X_eq**2)")
    
    # Substituting the calculated values into the equation
    lhs = eq_val
    a_val = R**(2.0/3.0)
    rhs = R * eq_val / (1 + a_val * eq_val**2)

    print("\nSubstituting the calculated values:")
    print(f"{lhs:.6f} = {R} * {eq_val:.6f} / (1 + {a_val:.6f} * {eq_val:.6f}**2)")
    print(f"{lhs:.6f} = {R * eq_val:.6f} / (1 + {a_val * eq_val**2:.6f})")
    print(f"{lhs:.6f} = {R * eq_val:.6f} / {1 + a_val * eq_val**2:.6f}")
    print(f"{lhs:.6f} ≈ {rhs:.6f}")


solve()
<<<1.049182>>>