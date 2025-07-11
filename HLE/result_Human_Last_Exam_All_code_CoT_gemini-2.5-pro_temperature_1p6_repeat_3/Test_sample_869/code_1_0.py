import math

def solve_probability():
    """
    Calculates the probability P = F/S based on the plan.
    S = 25! / (5!)^5
    F = 5!
    P = F / S
    """
    # Numerator of S
    s_numerator = math.factorial(25)

    # Denominator of S
    s_denominator = math.factorial(5)**5

    # Total number of distributions S
    # S = s_numerator / s_denominator

    # Number of favorable distributions F
    f_val = math.factorial(5)

    # The probability P is F / S, which simplifies to:
    # P = 5! / (25! / (5!)^5) = (5! * (5!)^5) / 25!
    # Let's calculate the numerator and denominator of this final fraction.
    
    p_numerator = f_val * s_denominator
    p_denominator = s_numerator

    print(f"The total number of items N = 25")
    print(f"The number of types T = 5, with 5 copies of each type.")
    print(f"The number of individuals is 5, each receiving 5 items.")
    print("-" * 30)
    print("Step 1: Calculate the total number of possible distributions (S)")
    print(f"S = 25! / (5! * 5! * 5! * 5! * 5!)")
    print(f"S = {s_numerator} / ({math.factorial(5)}^5)")
    print(f"S = {s_numerator} / {s_denominator}")
    S = s_numerator / s_denominator
    print(f"S = {S}")
    print("-" * 30)
    
    print("Step 2: Calculate the number of favorable distributions (F)")
    print("This is based on assigning a unique specialty type to each individual and assuming the only valid configuration is when each individual holds all 5 items of their specialty type.")
    print(f"F = 5! = {f_val}")
    print("-" * 30)
    
    print("Step 3: Calculate the probability P = F / S")
    print(f"P = F / S = (5!) / (25! / (5!)^5)")
    print("P = (5! * (5!)^5) / 25!")
    print(f"The final equation for the probability is:")
    print(f"P = ({f_val} * {math.factorial(5)}^5) / {s_numerator}")
    print(f"P = {p_numerator} / {p_denominator}")

    probability = p_numerator / p_denominator
    print(f"P ≈ {probability}")

solve_probability()