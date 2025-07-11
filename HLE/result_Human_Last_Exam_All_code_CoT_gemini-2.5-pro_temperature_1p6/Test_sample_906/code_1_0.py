import math

def solve_steady_state_pi0():
    """
    This function derives and prints the steady-state probability pi_0
    for the given birth-death process.
    """
    
    print("Derivation of the steady-state probability pi_0:")
    print("-" * 50)
    
    # Step 1: State the detailed balance equations for a birth-death process.
    print("Step 1: The detailed balance equations are pi_i * lambda_i = pi_{i+1} * mu_{i+1}.")
    print("         Given rates: lambda_i = lambda / (i+1), mu_{i+1} = mu.")
    print("         Substituting the rates: pi_i * (lambda / (i+1)) = pi_{i+1} * mu.")
    print("-" * 50)

    # Step 2: Express pi_k in terms of pi_0.
    print("Step 2: Rearrange to find pi_{i+1} in terms of pi_i.")
    print("         pi_{i+1} = pi_i * (lambda / mu) * (1 / (i+1)).")
    print("         Let rho = lambda / mu. The equation becomes: pi_{i+1} = pi_i * rho / (i+1).")
    print("\n         Let's find the general form for pi_k:")
    print("         For k=1: pi_1 = pi_0 * rho / 1 = pi_0 * rho")
    print("         For k=2: pi_2 = pi_1 * rho / 2 = (pi_0 * rho) * rho / 2 = pi_0 * rho^2 / 2!")
    print("         For k=3: pi_3 = pi_2 * rho / 3 = (pi_0 * rho^2 / 2!) * rho / 3 = pi_0 * rho^3 / 3!")
    print("\n         The general pattern is: pi_k = pi_0 * (rho^k / k!)")
    print("-" * 50)
    
    # Step 3: Use the normalization condition.
    print("Step 3: The sum of all probabilities must be 1: sum(pi_k for k=0 to inf) = 1.")
    print("         Substituting the formula for pi_k:")
    print("         sum(pi_0 * (rho^k / k!) for k=0 to inf) = 1")
    print("         pi_0 * sum(rho^k / k! for k=0 to inf) = 1")
    print("-" * 50)

    # Step 4: Recognize the series.
    print("Step 4: The sum is the Taylor series for the exponential function e^rho.")
    print("         sum(rho^k / k! for k=0 to inf) = e^rho")
    print("-" * 50)

    # Step 5: Solve for pi_0.
    print("Step 5: Substitute the series sum back into the equation:")
    print("         pi_0 * e^rho = 1")
    print("         Solving for pi_0 gives:")
    print("         pi_0 = 1 / e^rho = e^(-rho)")
    print("-" * 50)

    # Final Answer
    print("\nFinal Equation:")
    # Printing each "number" or component of the final equation e**(-rho)
    equation_parts = ["pi_0", "=", "e", "**", "(", "-", "rho", ")"]
    print(" ".join(equation_parts))

if __name__ == '__main__':
    solve_steady_state_pi0()