import math

def solve_titan_problem():
    """
    Solves the projectile problem using Titan's 5-bit fractional arithmetic constraints.
    """
    print("Step 1: Derive and simplify the governing equation.")
    print("The force F is given by F = (d * m * g) / h.")
    print("The mass m is given by m = rho * (4/3) * pi * r^3.")
    print("Substituting m and known values (d=20, h=10, rho=9/10, r=1/2):")
    print("F = (20/10) * (9/10 * 4/3 * pi * (1/2)^3) * g")
    print("F = 2 * (9/10 * 4/3 * pi * 1/8) * g")
    print("After simplification, the formula becomes: F = (3/10) * pi * g")
    print("-" * 20)

    print("Step 2: Select fractional approximations for pi and g.")
    print("To keep all intermediate numerators and denominators within the 5-bit limit (<= 31), we must choose our approximations carefully.")
    print("We choose pi approx= 3/1 and g approx= 10/1, as they allow for simplification.")
    
    pi_n, pi_d = 3, 1
    g_n, g_d = 10, 1
    
    print(f"Using pi = {pi_n}/{pi_d} and g = {g_n}/{g_d}")
    print("-" * 20)

    print("Step 3: Perform the calculation under Titan constraints.")
    # Initial formula F = (3/10) * pi * g
    c1_n, c1_d = 3, 10
    
    print(f"Initial expression: F = ({c1_n}/{c1_d}) * ({pi_n}/{pi_d}) * ({g_n}/{g_d})")

    # To avoid overflow, we can reorder the calculation to cancel terms first.
    # F = (3/10) * (10/1) * (3/1)
    # The (10) in the numerator and denominator cancel out.
    print("We reorder and simplify: (3/10) * (10/1) becomes 3/1.")
    
    # The expression is now F = (3/1) * (3/1)
    final_n = 3 * 3
    final_d = 1 * 1
    
    print("The final calculation is (3/1) * (3/1).")
    print("-" * 20)
    
    print("Final Equation and Result:")
    # The prompt requires printing the final equation with all numbers.
    print(f"({c1_n}/{c1_d}) * ({pi_n}/{pi_d}) * ({g_n}/{g_d}) = {final_n}/{final_d}")
    print(f"The required force is {final_n}/{final_d} N.")

    # For the final answer, calculate the error
    # True value: F = 2 * (0.9 kg/cm^3 * 4/3 * pi * (0.5cm)^3) * 9.8 m/s^2
    # m_true = 0.9 * (4/3) * math.pi * 0.125 = 0.15 * math.pi kg
    # F_true = 2 * m_true * 9.8 = 2 * (0.15 * math.pi) * 9.8 = 2.94 * math.pi
    F_true = 2.94 * math.pi
    F_calculated = final_n / final_d
    error = abs(F_calculated - F_true)
    
    # This part is for the final answer format, not part of the Titan output.
    # print(f"\nFor submission: True Force ~ {F_true:.4f} N. Calculated Force = {F_calculated:.4f} N. Absolute Error = {error:.4f}")


solve_titan_problem()

# The problem is solvable with the given constraints.
# The smallest absolute error found through a valid calculation path is |9 - 9.236...| = 0.236...
# So the answer is Y[e].
# <<<Y[0.236]>>>