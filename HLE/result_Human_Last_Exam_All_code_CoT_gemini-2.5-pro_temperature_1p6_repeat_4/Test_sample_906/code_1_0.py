import math

def solve_steady_state():
    """
    This function derives the steady-state probability pi_0 for the given birth-death process.
    """
    print("Step 1: The general formula for steady-state probabilities in a birth-death process is:")
    print("pi_k = pi_0 * product_{i=0 to k-1} (lambda_i / mu_{i+1})")
    print("-" * 60)

    print("Step 2: Substitute the given rates for this specific process:")
    print("Birth rate lambda_i = lambda / (i + 1)")
    print("Death rate mu_{i+1} = mu (for any i >= 0, mu_i is constant for i>0)")
    print("\nThis gives:")
    print("pi_k = pi_0 * product_{i=0 to k-1} ( (lambda / (i + 1)) / mu )")
    print("-" * 60)

    print("Step 3: Simplify the product term. Let rho = lambda / mu.")
    print("The product becomes: product_{i=0 to k-1} (rho / (i + 1))")
    print("Expanding the product: (rho/1) * (rho/2) * (rho/3) * ... * (rho/k)")
    print("This simplifies to: (rho^k) / (1 * 2 * 3 * ... * k) = rho^k / k!")
    print("\nSo, the relationship between pi_k and pi_0 is:")
    print("pi_k = pi_0 * (rho^k / k!)")
    print("-" * 60)

    print("Step 4: Apply the normalization condition that all probabilities must sum to 1:")
    print("sum_{k=0 to infinity} pi_k = 1")
    print("\nSubstitute the expression for pi_k:")
    print("sum_{k=0 to infinity} pi_0 * (rho^k / k!) = 1")
    print("\nFactor out pi_0:")
    print("pi_0 * (sum_{k=0 to infinity} rho^k / k!) = 1")
    print("-" * 60)

    print("Step 5: Recognize the summation as the Maclaurin (Taylor) series for e^x.")
    print("The series sum_{k=0 to infinity} x^k / k! is the expansion of e^x.")
    print("In our case, x = rho, so the sum is e^rho.")
    print("\nThe equation becomes:")
    print("pi_0 * e^rho = 1")
    print("-" * 60)
    
    print("Step 6: Solve for pi_0 to find the final expression.")
    print("Dividing both sides by e^rho, we get:")
    print("\nFinal Answer Equation: pi_0 = e^(-rho)")

solve_steady_state()