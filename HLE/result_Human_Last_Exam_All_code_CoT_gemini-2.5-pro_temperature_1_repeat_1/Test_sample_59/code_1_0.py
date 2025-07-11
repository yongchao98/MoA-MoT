import sys

def solve_probability_formula():
    """
    This function derives and prints the formula for the probability of a link
    in a jointly exchangeable random graph.
    """

    # Set encoding to handle unicode characters like the integral sign
    if sys.stdout.encoding != 'utf-8':
        sys.stdout.reconfigure(encoding='utf-8')

    print("The probability to draw a link y_ij in a jointly exchangeable random graph is derived from its representation via a random graphon.")
    print("The final probability is the expected value of the edge density of the random graphon.")
    print("\n--- FORMULA ---")

    # Define the components of the mathematical formula
    prob_term = "P(y_ij = 1)"
    expectation_term = "E_W"
    integral_sign = "\u222B"  # Unicode for the integral sign
    lower_bound = "0"
    upper_bound = "1"
    integrand = "W(u, v)"
    differentials = "du dv"

    # Construct the formula string for a more readable representation
    # P(y_ij = 1) = E_W[ ∫ (from 0 to 1) ∫ (from 0 to 1) W(u, v) du dv ]
    formula = (
        f"{prob_term} = {expectation_term}["
        f"{integral_sign} (from {lower_bound} to {upper_bound}) "
        f"{integral_sign} (from {lower_bound} to {upper_bound}) "
        f"{integrand} {differentials}]"
    )

    print(formula)

    print("\n--- EXPLANATION OF TERMS ---")
    print(f"1. P(y_ij = 1): The probability of a link existing between any two distinct nodes i and j.")
    print(f"2. W: The random graphon, which is a symmetric function W: [0,1]x[0,1] -> [0,1]. Its randomness is specified by the random measure F.")
    print(f"3. E_W: The expectation taken over the distribution of the random graphon W.")
    print(f"4. {integral_sign}: The integral sign. The double integral calculates the expected edge density for a specific realization of W.")
    print(f"   - The integration variables 'u' and 'v' represent the latent positions of the two nodes, drawn from the uniform measure U on [0,1].")
    print(f"\n--- THE NUMBERS IN THE EQUATION ---")
    print(f"The number '{lower_bound}' is the lower bound of integration for the latent node variables.")
    print(f"The number '{upper_bound}' is the upper bound of integration for the latent node variables.")
    print("These bounds define the latent space [0,1] where node attributes are generated according to the uniform measure U.")


solve_probability_formula()