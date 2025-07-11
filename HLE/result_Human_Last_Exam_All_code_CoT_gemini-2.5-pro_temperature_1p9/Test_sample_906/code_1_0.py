import math

def solve_steady_state_probability():
    """
    This function prints the step-by-step derivation for the steady-state
    probability pi_0 for the given birth-death process.
    """
    print("Derivation for the steady-state probability pi_0:\n")
    print("Let rho = lambda / mu\n")
    
    print("Step 1: The general formula for steady-state probability pi_n in a birth-death process is:")
    print("pi_n = pi_0 * PRODUCT_{i=0 to n-1} (lambda_i / mu_{i+1})\n")
    
    print("Step 2: Substitute the given rates into the formula.")
    print("Rates: lambda_i = lambda / (i+1) and mu_{i+1} = mu")
    print("The ratio is: (lambda / (i+1)) / mu = (lambda/mu) / (i+1) = rho / (i+1)\n")
    
    print("Step 3: Calculate the product term for pi_n.")
    print("pi_n = pi_0 * PRODUCT_{i=0 to n-1} (rho / (i+1))")
    print("pi_n = pi_0 * (rho/1) * (rho/2) * ... * (rho/n)")
    print("This simplifies to:")
    print("pi_n = pi_0 * (rho^n / n!)\n")
    
    print("Step 4: Use the normalization condition that the sum of all probabilities is 1.")
    print("SUM_{n=0 to infinity} pi_n = 1")
    print("SUM_{n=0 to infinity} pi_0 * (rho^n / n!) = 1")
    print("pi_0 * SUM_{n=0 to infinity} (rho^n / n!) = 1\n")
    
    print("Step 5: Recognize that the summation is the Taylor series for e^rho.")
    print("SUM_{n=0 to infinity} (rho^n / n!) = e^rho")
    print("The equation becomes: pi_0 * e^rho = 1\n")
    
    print("Step 6: Solve for pi_0.")
    print("The final expression for the steady-state probability pi_0 is:")
    print("pi_0 = 1 / e^rho = e^(-rho)")

# Execute the function to show the derivation.
solve_steady_state_probability()