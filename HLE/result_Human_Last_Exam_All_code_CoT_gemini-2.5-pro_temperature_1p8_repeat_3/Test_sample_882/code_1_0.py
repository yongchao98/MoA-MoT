import math

def solve_for_p():
    """
    Solves for the mixed strategy equilibrium probability p.
    """
    N = 8
    M = 3 * N

    print(f"Game parameters: N = {N} races, M = {M} players.")
    print("In a mixed strategy equilibrium, a player is indifferent between playing Discrete (D) and Continuous (C).")
    print("The equilibrium condition E_D(p) = E_C(p) leads to the equation for p:")
    print(f"p = (1 - (1 - p/{N})**{M}) / (1 - (p * ({N}-1)/{N})**{M})")
    print("We solve this equation using fixed-point iteration.")

    # Initial guess for p. We know from the problem description that p is likely close to 1.
    p = 0.98

    # Perform fixed-point iteration to find p
    for i in range(50):
        p_prev = p
        
        # Calculate terms of the equation g(p) = p
        numerator = 1 - (1 - p / N)**M
        denominator = 1 - (p * (N - 1) / N)**M
        
        if denominator == 0:
            # Avoid division by zero, although it's unlikely with p < 1
            p = float('inf')
        else:
            p = numerator / denominator
        
        # Check for convergence
        if abs(p - p_prev) < 1e-12:
            break
            
    print(f"\nSolved value for p: {p:.10f}")
    
    # Calculate the quantities required by the problem
    one_minus_p = 1 - p
    result_val = 10000 * one_minus_p
    final_answer = math.floor(result_val)
    
    print(f"Value of (1-p): {one_minus_p:.10f}")
    print(f"Value of 10000 * (1-p): {result_val:.10f}")
    print(f"The floor of 10000 * (1-p) is: {final_answer}")
    
    return final_answer

if __name__ == "__main__":
    solve_for_p()
