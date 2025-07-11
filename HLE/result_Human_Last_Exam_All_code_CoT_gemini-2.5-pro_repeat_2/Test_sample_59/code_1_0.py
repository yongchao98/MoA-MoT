import sympy

def calculate_link_probability():
    """
    Calculates the probability of a link in a jointly exchangeable random graph
    based on a specific latent variable model.
    """
    # Define the symbolic variable for our calculation
    f = sympy.symbols('f')

    # Explain the model and methodology
    print("To find the probability, we assume a specific generative model for a jointly exchangeable graph:")
    print("1. Each node `i` has a latent variable x_i ~ U[0,1].")
    print("2. A global random threshold F ~ U[0,1] is chosen.")
    print("3. A link y_ij exists if min(x_i, x_j) > F.")
    print("-" * 50)
    
    # Step 1: Conditional probability
    print("Step 1: Calculate the probability conditional on F=f.")
    print("P(y_ij=1 | F=f) = P(x_i > f and x_j > f)")
    print("                 = P(x_i > f) * P(x_j > f)   (due to independence)")
    print("                 = (1 - f) * (1 - f) = (1 - f)^2")
    print("-" * 50)

    # Step 2: Unconditional probability (Integration)
    print("Step 2: Average over all possible values of f by integrating from 0 to 1.")
    integrand = (1 - f)**2
    print(f"P(y_ij=1) = Integral from 0 to 1 of ({integrand}) df")
    
    # Perform the integration symbolically
    antiderivative = sympy.integrate(integrand, f)
    print(f"The antiderivative (indefinite integral) is: {antiderivative}")
    print("-" * 50)

    # Evaluate the definite integral and display the final equation components
    print("The final equation is the value of the antiderivative at the upper bound (f=1) minus its value at the lower bound (f=0).")
    
    # Numbers for the final equation
    upper_bound = 1
    lower_bound = 0
    upper_bound_val = antiderivative.subs(f, upper_bound)
    lower_bound_val = antiderivative.subs(f, lower_bound)
    final_probability = upper_bound_val - lower_bound_val
    
    print(f"\nEvaluating the final equation: [ {antiderivative} ] from {lower_bound} to {upper_bound}")
    print(f"Value at upper bound f={upper_bound}: {upper_bound_val}")
    print(f"Value at lower bound f={lower_bound}: {lower_bound_val}")
    print(f"Result = {upper_bound_val} - ({lower_bound_val})")
    
    print("\n" + "="*50)
    print(f"The final probability to draw a link y_ij is: {final_probability}")
    print("="*50)

if __name__ == "__main__":
    calculate_link_probability()
    # The final answer is 1/3
    # <<<1/3>>>
