def derive_pi0():
    """
    This script derives and prints the formula for the steady-state probability pi_0
    for the given birth-death process.
    """
    print("Step-by-step derivation of the steady-state probability pi_0:")
    print("=" * 60)
    
    # Introduction of terms
    print("The system is a birth-death process with the following transition rates:")
    print("  - Birth rate from state i to i+1: lambda_i = lambda / (i+1)")
    print("  - Death rate from state i to i-1: mu_i = mu")
    print("We define rho = lambda / mu.\n")

    # Step 1: Detailed balance equations
    print("Step 1: The detailed balance equations for a steady-state system are:")
    print("  pi_i * lambda_i = pi_{i+1} * mu_{i+1}")
    print("Substituting the given rates, we get:")
    print("  pi_i * (lambda / (i+1)) = pi_{i+1} * mu\n")

    # Step 2: Recurrence relation
    print("Step 2: Rearranging the equation to find a recurrence relation for pi_{i+1}:")
    print("  pi_{i+1} = pi_i * (lambda / (mu * (i+1)))")
    print("Using rho = lambda / mu, this simplifies to:")
    print("  pi_{i+1} = pi_i * (rho / (i+1))\n")

    # Step 3: Express pi_i in terms of pi_0
    print("Step 3: Expressing each pi_i in terms of pi_0 by induction:")
    print("  For i=0: pi_1 = pi_0 * (rho / 1) = pi_0 * rho")
    print("  For i=1: pi_2 = pi_1 * (rho / 2) = (pi_0 * rho) * (rho / 2) = pi_0 * (rho^2 / 2!)")
    print("  For i=2: pi_3 = pi_2 * (rho / 3) = (pi_0 * rho^2 / 2!) * (rho / 3) = pi_0 * (rho^3 / 3!)")
    print("The general formula is therefore:")
    print("  pi_i = pi_0 * (rho^i / i!)\n")

    # Step 4: Normalization condition
    print("Step 4: Using the normalization condition that the sum of all probabilities is 1:")
    print("  sum_{i=0 to infinity} pi_i = 1")
    print("Substituting our formula for pi_i:")
    print("  sum_{i=0 to infinity} (pi_0 * (rho^i / i!)) = 1")
    print("Factoring out pi_0:")
    print("  pi_0 * (sum_{i=0 to infinity} rho^i / i!) = 1\n")
    
    # Step 5: Recognize the Taylor series
    print("Step 5: The sum is the well-known Taylor series expansion for the exponential function e^x, where x = rho:")
    print("  sum_{i=0 to infinity} rho^i / i! = e^rho")
    print("So, the equation becomes:")
    print("  pi_0 * e^rho = 1\n")

    # Step 6: Solve for pi_0
    print("Step 6: Solving for pi_0 gives the final expression:")
    print("  pi_0 = 1 / e^rho\n")

    # Final Answer
    print("The final equation for the steady-state probability pi_0 is:")
    print("  pi_0 = e^(-rho)")

    # Outputting numbers in the final equation as per the instructions
    print("\nThe numbers present in the final symbolic equation 'pi_0 = e^(-rho)' are:")
    print("The subscript of pi: 0")
    print("The implicit coefficient of rho in the exponent: -1")


if __name__ == '__main__':
    derive_pi0()
