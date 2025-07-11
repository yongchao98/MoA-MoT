import math

def solve_steady_state():
    """
    This function derives the steady-state probability pi_0 for the given
    birth-death process and prints the derivation steps.
    """
    print("Step 1: Express pi_n in terms of pi_0 using the detailed balance equations.")
    print("The birth rates are lambda_i = lambda / (i + 1).")
    print("The death rates are mu_i = mu.")
    print("\nThe balance equations are: pi_{i+1} * mu_{i+1} = pi_i * lambda_i")
    print("pi_{i+1} = pi_i * (lambda_i / mu_{i+1})")
    print("Given lambda_i = lambda / (i + 1) and mu_{i+1} = mu, we have:")
    print("pi_{i+1} = pi_i * ( (lambda / (i + 1)) / mu ) = pi_i * (lambda / ((i + 1) * mu))")
    print("\nLet's find the general form for pi_n in terms of pi_0:")
    print("pi_1 = pi_0 * (lambda / (1 * mu))")
    print("pi_2 = pi_1 * (lambda / (2 * mu)) = pi_0 * (lambda / (1 * mu)) * (lambda / (2 * mu)) = pi_0 * (lambda^2 / (2 * 1 * mu^2)) = pi_0 * (lambda^2 / (2! * mu^2))")
    print("pi_3 = pi_2 * (lambda / (3 * mu)) = pi_0 * (lambda^2 / (2! * mu^2)) * (lambda / (3 * mu)) = pi_0 * (lambda^3 / (3! * mu^3))")
    print("\nBy induction, we can see the pattern is:")
    print("pi_n = pi_0 * (lambda^n / (n! * mu^n))")
    print("Let rho = lambda / mu. The expression becomes:")
    print("pi_n = pi_0 * (rho^n / n!)")
    
    print("\n" + "="*50 + "\n")
    
    print("Step 2: Use the normalization condition Sum(pi_n for n=0 to inf) = 1.")
    print("Sum(pi_0 * (rho^n / n!)) = 1")
    print("pi_0 * Sum(rho^n / n!) = 1")
    
    print("\n" + "="*50 + "\n")

    print("Step 3: Recognize the series and solve for pi_0.")
    print("The series Sum(rho^n / n!) for n from 0 to infinity is the Taylor series for e^rho.")
    print("So, the equation becomes:")
    print("pi_0 * e^rho = 1")
    
    print("\n" + "="*50 + "\n")

    print("Step 4: Final expression for pi_0.")
    print("Solving for pi_0, we get:")
    print("pi_0 = 1 / e^rho")
    print("The final equation is:")
    print("pi_0 = e^(-rho)")

solve_steady_state()