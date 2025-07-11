import math

def calculate_limiting_probability():
    """
    Calculates the limiting probability P(n) as n approaches infinity.

    The problem states that replication is successful if for all 8 bases (x=0 to 7),
    the count of that base (c_x) does not satisfy c_x mod 8 = x.

    As n -> infinity, the distribution of c_x mod 8 becomes uniform over {0, ..., 7}.
    Therefore, the probability that c_x mod 8 = x approaches 1/8 for any given x.
    The probability of the success condition for one base (c_x mod 8 != x) is 1 - 1/8 = 7/8.

    Since the events for each base are independent in the limit, the total probability is
    the product of the individual probabilities for all 8 bases.
    """

    # The numbers that constitute the final equation for the probability
    base_numerator = 7
    base_denominator = 8
    exponent = 8

    # Calculate the final numerator and denominator from the closed-form expression
    final_numerator = base_numerator**exponent
    final_denominator = base_denominator**exponent

    # Calculate the final probability as a decimal value
    probability = final_numerator / final_denominator

    # Output the components of the equation and the final result
    print("The closed-form expression for the limiting probability is (7/8)^8.")
    print("\n--- Equation Breakdown ---")
    print(f"Base of the fraction: {base_numerator}/{base_denominator}")
    print(f"Exponent: {exponent}")
    print(f"Final equation: P = ({base_numerator}/{base_denominator})^{exponent} = {final_numerator} / {final_denominator}")
    
    print("\n--- Final Answer ---")
    print(f"The limiting probability as a fraction is: {final_numerator}/{final_denominator}")
    print(f"The limiting probability as a decimal is: {probability}")

# Execute the function to find and print the solution
calculate_limiting_probability()