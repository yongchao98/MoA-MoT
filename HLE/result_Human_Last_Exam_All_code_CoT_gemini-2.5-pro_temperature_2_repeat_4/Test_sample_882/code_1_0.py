import math

def solve_equilibrium_probability():
    """
    This function solves for the value of x = 1-p based on the game's 
    Nash equilibrium condition, and then computes the final requested value.

    The game has N races and M=3N players. For N=8, M=24.
    The problem states that for N=8, a mixed strategy Nash Equilibrium exists
    where players play the 'discrete strategy' with probability p.

    The expected payoff for any player in this symmetric game must be N/M = 1/3.
    In a mixed-strategy equilibrium, the expected payoff for any strategy played
    with a positive probability must be equal to this value.

    The expected payoff for the discrete strategy is E_D = (1 - (1 - p/N)^M) / (M*p/N).
    Setting E_D = N/M = 1/3, we get the equation for p:
      (1 - (1 - p/8)^24) / (24*p/8) = 1/3
    This simplifies to:
      1 - (1 - p/8)^24 = p

    To find the required value, we solve for x = 1-p. Substituting p = 1-x:
      1 - (1 - (1-x)/8)^24 = 1-x
      x = (1 - (1-x)/8)^24
      x = ((8 - 1 + x)/8)^24
      x = ((7+x)/8)^24
      
    We will solve this equation for x numerically using fixed-point iteration.
    """
    
    # Constants from the problem
    N = 8
    M = 3 * N
    
    # The final equation to solve for x = 1 - p
    # is: x = ((7 + x) / 8) ^ 24
    
    print("The problem parameters are:")
    print(f"N (number of races) = {N}")
    print(f"M (number of players) = {M}")
    
    print("\nThe equation for p is derived from setting the expected payoff of the discrete strategy to 1/3:")
    print(f"1 - (1 - p / {N})^{M} = p")
    print("Which for N=8 and M=24 becomes:")
    print(f"1 - (1 - p / {N})^{{{M}}} = p")
    print("\nSubstituting x = 1 - p, the equation to solve for x becomes:")
    num1, num2, num3 = 7, 8, 24
    print(f"x = (({num1} + x) / {num2})^{{{num3}}}")

    # Initial guess for x. From theoretical analysis, we know x is a small positive number.
    x = 0.04
    
    # Perform fixed-point iteration until convergence
    for _ in range(100):
        x_new = ((num1 + x) / num2)**num3
        if abs(x_new - x) < 1e-15:
            x = x_new
            break
        x = x_new
        
    print(f"\nSolved value for x (which is 1-p) is: {x:.10f}")

    # Calculate the final result
    result = math.floor(10000 * x)
    
    print(f"\nThe expression to calculate is floor(10000 * (1-p)) = floor(10000 * x).")
    print(f"10000 * x = {10000 * x:.6f}")
    print(f"The final answer is: {result}")
    
    return result

if __name__ == "__main__":
    solve_equilibrium_probability()