import numpy as np

def compute_integral_value():
    """
    This function computes the value of the integral based on the derived symbolic expression:
    Value = (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
    """
    pi = np.pi

    # Calculate each term of the final expression
    term1 = (8 / 15) * (pi**8)
    term2 = (1 / 3) * (pi**2)
    term3 = - (1 / 2) * pi
    term4 = 1.0

    # Sum the terms to get the final value
    total_value = term1 + term2 + term3 + term4

    # Print the final equation and the value of each number/term in it
    print("The final simplified expression for the integral is:")
    print("I = (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1")
    print("\nEvaluating each term:")
    print(f"The number for the pi^8 term is 8/15, and the term's value is: {term1}")
    print(f"The number for the pi^2 term is 1/3, and the term's value is: {term2}")
    print(f"The number for the pi term is -1/2, and the term's value is: {term3}")
    print(f"The constant term is 1, and its value is: {term4}")
    
    # Print the final numerical result
    print("\n------------------------------------")
    print(f"The total numerical value is: {total_value}")
    print("------------------------------------")

# Execute the function
compute_integral_value()