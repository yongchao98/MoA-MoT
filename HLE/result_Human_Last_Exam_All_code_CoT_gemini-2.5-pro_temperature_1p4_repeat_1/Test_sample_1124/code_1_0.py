import math

def main():
    """
    This script defines a modified logistic map and verifies its equilibrium point
    for a specific value of R.
    """
    # The value of R from the problem description
    R = 3.57

    # Define the numerator and denominator for the constant in our modified map
    numerator = 4
    denominator = 3
    constant = float(numerator) / denominator

    print("This script demonstrates a modified logistic map with a specific equilibrium point.")
    print("-" * 70)

    # Output the equation and its components as requested
    print(f"The standard logistic map is X_n+1 = R * X_n * (1 - X_n).")
    print(f"It is modified by changing the '1' to '{numerator}/{denominator}'.")
    print(f"The new equation is: X_n+1 = R * X_n * ({numerator}/{denominator} - X_n)")
    print(f"We will analyze this map for the given parameter R = {R}")
    print("-" * 70)

    # Step 1: Analytically calculate the equilibrium point X*
    # The fixed-point equation is X* = R * X* * (4/3 - X*). For X* != 0:
    # 1 = R * (4/3 - X*), which gives X* = 4/3 - 1/R.
    x_equilibrium = constant - (1.0 / R)

    print("Step 1: Calculate the analytical equilibrium point (X*).")
    print(f"The formula for the equilibrium point is: X* = {numerator}/{denominator} - 1/R")
    print(f"Plugging in R = {R}: X* = {constant} - {1.0 / R}")
    print(f"The calculated equilibrium point is: X* = {x_equilibrium}")
    print("\nThis value is approximately 1.05, as requested by the prompt.")
    print("-" * 70)

    # Step 2: Verify that f(X*, R) = X* by plugging X* back into the map
    # An equilibrium point X* must satisfy the condition: f(X*) = X*
    
    # Define the map function for clarity
    def modified_logistic_map(x, r):
        return r * x * (constant - x)

    # Apply the map to the calculated equilibrium point
    x_next = modified_logistic_map(x_equilibrium, R)

    print("Step 2: Verify that this is a true equilibrium point.")
    print("We must show that if we input X_n = X*, the output X_n+1 is also X*.")
    print("The equation with all numbers plugged in is:")
    print(f"X_n+1 = {R} * {x_equilibrium} * ({constant} - {x_equilibrium})")
    print(f"The result of this calculation is: X_n+1 = {x_next}")
    print("-" * 70)

    # Final conclusion
    print("Conclusion:")
    print(f"Calculated equilibrium point X* = {x_equilibrium}")
    print(f"Result of map f(X*, R)     = {x_next}")
    
    # Check if the input and output are equal (using a tolerance for float precision)
    if math.isclose(x_equilibrium, x_next):
        print("\nSince the output is equal to the input, the point is a confirmed equilibrium point.")
    else:
        print("\nThe output is not equal to the input. There is an error.")

if __name__ == "__main__":
    main()