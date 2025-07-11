import math

def solve_alpha():
    """
    Solves for the value of alpha based on the problem description.
    
    The plan is as follows:
    1.  The solvability condition for the perturbed system is x_0^4 + y_0^4 = C,
        where C = 3*alpha / (1 - exp(-3*T)).
    2.  The integral is I = integral from 0 to A of (C - x_0^4)*x_0^3 dx.
    3.  Evaluating the integral gives I = (A^4/8)*(2*C - A^4).
    4.  The problem statement defines A and C in a way that leads to an inconsistent, unsolvable
        polynomial for alpha. The most reasonable interpretation is that the integration limit A
        is the natural boundary of the problem, which is C^(1/4).
    5.  Assuming A^4 = C, the integral simplifies to I = C^2 / 8.
    6.  We set I = B, which is given, to find the value of C.
    7.  From the value of C, we calculate alpha.
    """
    
    # Given constants
    T = math.log(10)
    B = 0.5 * (10**20) / (99**2)
    
    # Step 1: Solve for C from I = C^2 / 8 = B
    # C^2 = 8 * B
    C_squared = 8 * B
    C = math.sqrt(C_squared)
    
    # Step 2: Solve for alpha from C = 3*alpha / (1 - exp(-3*T))
    # alpha = C * (1 - exp(-3*T)) / 3
    term_T = 1 - math.exp(-3 * T)
    alpha = C * term_T / 3
    
    # Output the results
    print("Based on the interpretation that A^4 = C, we solve C^2 / 8 = B.")
    print(f"Given T = ln(10) = {T:.4f}")
    print(f"Given B = 1/2 * 10^20 / 99^2 = {B:.4e}")
    print(f"Calculated C = sqrt(8 * B) = {C:.4e}")
    
    print("\nThen we solve for alpha using C = 3 * alpha / (1 - exp(-3*T)).")
    print("The final equation is:")
    print(f"{C:.4e} = (3 * alpha) / (1 - e^(-3 * {T:.4f}))")
    print(f"{C:.4e} = (3 * alpha) / {term_T}")
    print(f"alpha = {C:.4e} * {term_T} / 3")
    
    print(f"\nThe calculated value of alpha is: {alpha}")
    print(f"As a fraction, this is (74/11) * 10^7.")

solve_alpha()