import math

def modified_logistic_map():
    """
    This function implements a modified logistic map that creates a stable
    equilibrium point at X = 1.05 for R = 3.57.
    """
    # Parameter R, at which the standard map is chaotic
    R = 3.57
    
    # The desired equilibrium point, which is now part of the equation's structure
    X_eq = 1.05
    
    # An initial value for X, can be anything reasonable
    X_current = 0.5
    
    # Number of iterations to show convergence
    iterations = 50

    print("Running the modified logistic map simulation...")
    print(f"Goal: Achieve a stable equilibrium at X = {X_eq} with R = {R}.")
    print(f"Initial X value: {X_current}")
    print("-" * 50)

    # Iterate the map to show convergence
    for i in range(iterations):
        # Our new map equation: X_n+1 = X_eq - R * (X_n - X_eq)^2
        X_next = X_eq - R * (X_current - X_eq)**2
        
        # We don't need to print every single iteration, just enough to see the trend.
        if i < 10 or i == iterations - 1:
             print(f"Iteration {i+1:2d}: X = {X_next:.8f}")
        elif i == 10:
             print("...")

        # Stop if it gets too close to prevent floating point noise
        if math.isclose(X_current, X_next):
            print(f"\nConvergence reached at iteration {i+1}.")
            X_current = X_next
            break
            
        X_current = X_next
        
    print("-" * 50)
    print(f"The final value of X is approximately: {X_current:.8f}")
    
    print("\nThe modified logistic map equation that achieves this is:")
    # Print each number in the final equation as requested
    print(f"X_n+1 = {X_eq} - {R} * (X_n - {X_eq})^2")

# Run the simulation
modified_logistic_map()
