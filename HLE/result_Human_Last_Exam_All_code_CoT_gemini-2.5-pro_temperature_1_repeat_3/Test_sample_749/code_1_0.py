import math

def calculate_critical_h():
    """
    Calculates the critical branching probability h for a specific
    2-site reactor configuration to become supercritical.
    """
    
    # The equation for the critical threshold is h > n / (d1 * (d2 + sqrt(d3))) - s
    n = 5.0
    d1 = 2.0
    d2 = 1.0
    d3 = 2.0
    s = 1.0
    
    sqrt_d3 = math.sqrt(d3)
    
    denominator = d1 * (d2 + sqrt_d3)
    
    critical_h = n / denominator - s
    
    print("To solve for the critical value of h, we analyze the condition for a 'reactor' to be self-sustaining.")
    print("One of the most favorable reactor configurations involves sites -1 and 0 being red.")
    print("This leads to a 2-type branching process with a mean matrix M = (1+h) * M_0.")
    print("The system is supercritical if the largest eigenvalue of M is > 1.")
    print("This gives the condition: (1+h) * (2/5)*(1 + sqrt(2)) > 1")
    print("Rearranging for h, we get the inequality: h > 5 / (2 * (1 + sqrt(2))) - 1")
    print("\nCalculating the values from the final equation:")
    print(f"Numerator: {n}")
    print(f"Denominator term 1: {d1}")
    print(f"Denominator term 2: {d2}")
    print(f"Denominator term 3 (under square root): {d3}")
    print(f"Value of sqrt({d3}): {sqrt_d3}")
    print(f"Subtrahend: {s}")
    
    print(f"\nThe critical value h_c is: {critical_h}")
    
    print("\nSince the critical value h_c is a positive number, for any h -> 0, this reactor (and any other) will be subcritical.")
    print("Therefore, no self-sustaining colony can form near site 0.")
    print("The probability of infinitely many particles visiting site 0 must therefore tend to 0.")

calculate_critical_h()