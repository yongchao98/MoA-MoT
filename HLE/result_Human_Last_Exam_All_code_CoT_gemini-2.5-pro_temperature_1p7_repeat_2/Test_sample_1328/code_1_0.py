import math

def find_optimal_b(P, P1, P2, alpha):
    """
    Calculates the optimal feedback adjustment factor 'b'.

    Args:
        P (float): Half of the total power budget.
        P1 (float): Power used in the first transmission.
        P2 (float): Power used in the second transmission signal component.
        alpha (float): Weather-induced correlation.
    """
    print(f"--- Calculating Optimal b for P={P}, P1={P1}, P2={P2}, alpha={alpha} ---")

    # Step 1: Check if the power allocation is feasible
    power_remainder = 2 * P - P1 - P2
    if power_remainder < 0:
        print("Error: Power constraint cannot be met because P1 + P2 > 2P.")
        print(f"       {P1} + {P2} = {P1 + P2}, which is greater than 2*P = {2*P}.")
        return

    # Step 2: Calculate the boundary of the feasible region for b
    B = math.sqrt(power_remainder)

    # Step 3: Determine the optimal b based on the sign of alpha
    if alpha > 0:
        optimal_b = -B
        print("Analysis: Since alpha > 0, the objective function is maximized at the negative boundary of the feasible region for b.")
        print("The formula for optimal b is: b_opt = -sqrt(2*P - P1 - P2)")
        print(f"Calculation: b_opt = -sqrt(2*{P} - {P1} - {P2})")
        print(f"             b_opt = -sqrt({power_remainder})")
        print(f"Final Answer: b_opt = {optimal_b}")
    else:  # This branch handles both alpha < 0 and alpha = 0
        optimal_b = B
        if alpha < 0:
            print("Analysis: Since alpha < 0, the objective function is maximized at the positive boundary of the feasible region for b.")
        else:  # alpha == 0
            print("Analysis: Since alpha = 0, both boundaries are optimal. We choose the positive one by convention.")
        
        print("The formula for optimal b is: b_opt = sqrt(2*P - P1 - P2)")
        print(f"Calculation: b_opt = sqrt(2*{P} - {P1} - {P2})")
        print(f"             b_opt = sqrt({power_remainder})")
        print(f"Final Answer: b_opt = {optimal_b}")

# --- Example Usage ---
# Case 1: alpha is positive
find_optimal_b(P=20, P1=10, P2=5, alpha=0.6)
print("\n" + "="*70 + "\n")

# Case 2: alpha is negative
find_optimal_b(P=20, P1=10, P2=5, alpha=-0.3)
print("\n" + "="*70 + "\n")

# Case 3: alpha is zero
find_optimal_b(P=20, P1=10, P2=5, alpha=0.0)
print("\n" + "="*70 + "\n")

# Case 4: Infeasible power allocation
find_optimal_b(P=20, P1=15, P2=30, alpha=0.5)
