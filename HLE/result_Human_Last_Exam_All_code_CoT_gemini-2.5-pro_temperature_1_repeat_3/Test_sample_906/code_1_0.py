import math

def solve_and_explain():
    """
    This function explains the derivation of the steady-state probability pi_0
    for the given birth-death process and prints the final result.
    """
    
    print("Step 1: Define the detailed balance equations.")
    print("For a birth-death process in steady state, the rate of transitioning from state i-1 to i must equal the rate of transitioning from i to i-1.")
    print("This gives us the equation: pi_{i-1} * q_{i-1, i} = pi_i * q_{i, i-1} for i >= 1.")
    print("\nGiven the transition rates:")
    print("  - Birth rate (i -> i+1): q_{i, i+1} = lambda / (i + 1)")
    print("  - Death rate (i -> i-1): q_{i, i-1} = mu")
    print("\nThe balance equation is:")
    print("pi_{i-1} * (lambda / i) = pi_i * mu")

    print("\nStep 2: Find a recurrence relation for pi_i.")
    print("We rearrange the balance equation to solve for pi_i in terms of pi_{i-1}:")
    print("pi_i = (lambda / (i * mu)) * pi_{i-1}")
    print("Using the definition rho = lambda / mu, we get:")
    print("pi_i = (rho / i) * pi_{i-1}")

    print("\nStep 3: Express pi_i in terms of pi_0.")
    print("We can solve this recurrence relation by iterating from i=1:")
    print("  - For i = 1: pi_1 = (rho / 1) * pi_0")
    print("  - For i = 2: pi_2 = (rho / 2) * pi_1 = (rho / 2) * (rho * pi_0) = (rho**2 / 2!) * pi_0")
    print("  - For i = 3: pi_3 = (rho / 3) * pi_2 = (rho / 3) * (rho**2 / 2!) * pi_0 = (rho**3 / 3!) * pi_0")
    print("The general pattern is:")
    print("pi_i = (rho**i / i!) * pi_0")

    print("\nStep 4: Use the normalization condition to solve for pi_0.")
    print("The sum of all probabilities must be 1: Sum_{i=0 to infinity} pi_i = 1.")
    print("Substituting our expression for pi_i:")
    print("Sum_{i=0 to infinity} (rho**i / i!) * pi_0 = 1")
    print("pi_0 * (Sum_{i=0 to infinity} rho**i / i!) = 1")

    print("\nStep 5: Identify the infinite series.")
    print("The sum is the Taylor series expansion for the exponential function e^x, where x = rho:")
    print("Sum_{i=0 to infinity} rho**i / i! = e**rho")

    print("\nStep 6: Find the final expression for pi_0.")
    print("Substituting the sum back into the equation gives:")
    print("pi_0 * e**rho = 1")
    print("Solving for pi_0, we get:")
    print("\npi_0 = 1 / e**rho = e**(-rho)")
    
solve_and_explain()