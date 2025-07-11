import math

def solve_for_p():
    """
    Solves the equation p = 1 - (1 - p/N)^(3N) for N=8 numerically.
    """
    N = 8
    three_N = 3 * N

    # The equation for the equilibrium probability p is:
    # p = 1 - (1 - p/N)^(3N)
    # We can solve this using fixed-point iteration.
    # p_{i+1} = 1 - (1 - p_i/N)^(3N)
    
    # Start with an initial guess for p (must be between 0 and 1)
    # p=0 is a trivial solution, we seek the non-zero one.
    p = 0.5 
    
    # Iterate to find the fixed point, which is the value of p.
    # 30 iterations are more than enough for convergence to high precision.
    for i in range(30):
        p = 1 - (1 - p / N) ** three_N

    return p

def main():
    """
    Calculates the final value based on the equilibrium probability p.
    """
    N = 8
    three_N = 3 * N
    
    print(f"Finding the mixed strategy equilibrium for N = {N} players.")
    print(f"The equilibrium probability p is the non-zero solution to the equation:")
    print(f"p = 1 - (1 - p/{N})**{three_N}")
    print("-" * 30)

    p = solve_for_p()
    
    print(f"Numerically solved value for p: {p:.6f}")
    
    one_minus_p = 1 - p
    print(f"The value of (1-p) is: {one_minus_p:.6f}")

    result = 10000 * one_minus_p
    final_answer = math.floor(result)
    
    print("-" * 30)
    print("Calculating the final result as per the problem statement:")
    print(f"The expression is: floor(10000 * (1-p))")
    print(f"Substituting the value of p: floor(10000 * {one_minus_p:.6f})")
    print(f"This evaluates to: floor({result:.4f})")
    print(f"The final answer is: {final_answer}")

if __name__ == "__main__":
    main()