import math

def solve_steady_state():
    """
    This script derives the steady-state probability pi_0 for the given birth-death process.
    It prints the derivation steps and the final formula.
    """

    print("Step 1: The Birth-Death Process Formulation")
    print("------------------------------------------")
    print("The transition rates for the process are:")
    print("Birth rate (from state i to i+1):   lambda_i = lambda / (i+1)")
    print("Death rate (from state i to i-1):   mu_i = mu")
    print("Let the traffic intensity be defined as rho = lambda / mu.")
    print("\nThe detailed balance equations for the steady-state probabilities pi_i are:")
    print("pi_i * lambda_i = pi_{i+1} * mu_{i+1}")
    print("Rearranging for pi_{i+1}: pi_{i+1} = pi_i * (lambda_i / mu_{i+1})\n")

    print("Step 2: Expressing pi_k in terms of pi_0")
    print("------------------------------------------")
    print("We can find a general expression for pi_k in terms of pi_0 by recursive substitution:")
    print("pi_1 = pi_0 * (lambda_0 / mu_1) = pi_0 * ( (lambda/(0+1)) / mu ) = pi_0 * (rho / 1!)")
    print("pi_2 = pi_1 * (lambda_1 / mu_2) = [pi_0 * rho] * ( (lambda/(1+1)) / mu ) = pi_0 * (rho**2 / 2!)")
    print("pi_3 = pi_2 * (lambda_2 / mu_3) = [pi_0 * rho**2 / 2!] * ( (lambda/(2+1)) / mu ) = pi_0 * (rho**3 / 3!)")
    print("\nBy observing the pattern, we can write the general formula:")
    print("pi_k = pi_0 * (rho**k / k!)\n")

    print("Step 3: Applying the Normalization Condition")
    print("---------------------------------------------")
    print("The sum of all steady-state probabilities must be equal to 1:")
    print("Sum(pi_k for k=0 to infinity) = 1")
    print("\nSubstituting our expression for pi_k:")
    print("Sum(pi_0 * (rho**k / k!) for k=0 to infinity) = 1")
    print("pi_0 * Sum(rho**k / k! for k=0 to infinity) = 1\n")

    print("Step 4: Solving for pi_0")
    print("--------------------------")
    print("The sum is the well-known Taylor series expansion for the exponential function e^x, where x = rho:")
    print("Sum(rho**k / k! for k=0 to infinity) = e**rho")
    print("\nSubstituting this back into our normalization equation:")
    print("pi_0 * e**rho = 1")
    print("\nTherefore, the final expression for pi_0 is found by solving for it:")
    
    final_equation_expression = "pi_0 = e**(-rho)"
    print(f"Final Equation: {final_equation_expression}")

    print("\nIn this final equation:")
    print("pi_0 is the steady-state probability of being in state 0.")
    print("e is Euler's number, the base of the natural logarithm, approximately {:.4f}.".format(math.e))
    print("The number in the exponent of e is: -1")
    print("rho is the ratio lambda/mu.")

solve_steady_state()