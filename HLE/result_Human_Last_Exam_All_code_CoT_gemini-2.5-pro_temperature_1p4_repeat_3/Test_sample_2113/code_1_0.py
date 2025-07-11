import math

def solve_and_print():
    """
    Solves the problem by deriving the analytical forms of the trajectory and
    the entanglement echo, then computes the final value.
    """
    # Introduction to the method
    print("This problem is solved analytically based on the properties of the system.")
    print("1. The particle trajectory in a harmonic potential follows classical motion: z(t) = z(0) * cos(t).")
    print("2. The entanglement echo y(t) is solved from the integral equation, yielding y(t) = sin(2t) / (cos(2t))^(3/2).")
    print("\nWe will now calculate the value of (z1(t) / y(t))^2 at t = pi/8.")

    # Step 1: Define the specific time t
    t_val = math.pi / 8

    # Step 2: Define and calculate z1 at t
    # The trajectory z1(t) has initial condition z1(0) = 1.
    # Therefore, z1(t) = cos(t).
    z1_val = math.cos(t_val)

    # Step 3: Define and calculate y at t
    # y(t) = sin(2t) / (cos(2t))^(3/2)
    sin_2t = math.sin(2 * t_val) # sin(pi/4)
    cos_2t = math.cos(2 * t_val) # cos(pi/4)
    y_val = sin_2t / (cos_2t**1.5)

    # Step 4: Compute the final result
    result = (z1_val / y_val)**2
    
    # --- Output the numerical components of the final equation ---
    print("\n--- Calculation Steps ---")
    print(f"The value of t is pi/8.")
    
    print(f"\n1. Calculate z1(pi/8):")
    print(f"z1(pi/8) = cos(pi/8) = {z1_val:.6f}")

    print(f"\n2. Calculate y(pi/8):")
    print(f"y(pi/8) = sin(pi/4) / (cos(pi/4))^(3/2) = {y_val:.6f}")

    print("\n3. Compute the final expression: (z1(pi/8) / y(pi/8))^2")
    
    # --- Final Equation with numbers ---
    print("\nFinal Equation:")
    print(f"({z1_val:.4f} / {y_val:.4f})^2 = {result:.6f}")

    # For context, provide the exact analytical result
    analytical_result_str = "(sqrt(2) + 1) / 4"
    print(f"The exact analytical result is {analytical_result_str}.")
    
    return result

# Execute the solution
final_answer = solve_and_print()
# print(f'<<<{final_answer:.4f}>>>')
# Let's provide the final answer as a number outside the code block
# (sqrt(2) + 1)/4 = (1.41421356 + 1)/4 = 0.60355339
# Using 4 decimal places gives 0.6036

if __name__ == '__main__':
    # The script is designed to be self-explanatory when run.
    pass
