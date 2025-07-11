import math
from fractions import Fraction

def solve_race_growth_rate():
    """
    Calculates the difference between optimal and actual Kelly growth rates
    for a three-competitor race with mistaken probabilities.
    """
    # Step 1: Define problem parameters
    p_actual = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 4)]
    b = [4, 3, 3]  # Payout ratios b:1
    p_believed = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)]

    # Step 2: Calculate the optimal strategy (f_star) and optimal growth rate (W_star)
    # The edge for B and C is (1/4)*3 - (3/4) = 0, so f_B* = f_C* = 0.
    # The edge for A is (1/2)*4 - 1/2 = 1.5 > 0.
    # f_A* = edge / odds = 1.5 / 4 = 3/8.
    f_star = [Fraction(3, 8), Fraction(0), Fraction(0)]

    # Calculate wealth factors for W_star based on f_star
    r_a_star = 1 + b[0]*f_star[0] - f_star[1] - f_star[2]
    r_b_star = 1 - f_star[0] + b[1]*f_star[1] - f_star[2]
    r_c_star = 1 - f_star[0] - f_star[1] + b[2]*f_star[2]

    # Calculate W_star using actual probabilities
    # W* = (1/2)log(5/2) + (1/4)log(5/8) + (1/4)log(5/8) = log(5/4)
    w_star_val = math.log(Fraction(5,4))
    
    # Step 3: Calculate the mistaken betting fractions (f_actual)
    # These fractions are derived by maximizing the growth rate using the *believed*
    # probabilities. Solving the system of equations yields f_A=7/44, f_B=17/44.
    f_actual = [Fraction(7, 44), Fraction(17, 44), Fraction(0)]

    # Step 4: Calculate the actual growth rate (W) using f_actual and p_actual
    # Calculate wealth factors for W based on f_actual
    r_a_actual = 1 + b[0]*f_actual[0] - f_actual[1] - f_actual[2]
    r_b_actual = 1 - f_actual[0] + b[1]*f_actual[1] - f_actual[2]
    r_c_actual = 1 - f_actual[0] - f_actual[1] + b[2]*f_actual[2]

    # Calculate W using actual probabilities
    w_actual_val = p_actual[0]*math.log(r_a_actual) + p_actual[1]*math.log(r_b_actual) + p_actual[2]*math.log(r_c_actual)

    # Step 5: Print the derivation and final results
    print("Derivation of the difference W* - W:")
    print(f"W* = log({r_a_star / r_b_star}) = log(5/4)")
    print(f"W = {p_actual[0]}*log({r_a_actual}) + {p_actual[1]}*log({r_b_actual}) + {p_actual[2]}*log({r_c_actual})")
    print("W* - W = log(5/4) - [ (1/2)log(5/4) + (1/4)log(2) + (1/4)log(5/11) ]")
    print("         = (1/2)log(5/4) - (1/4)log(2) - (1/4)log(5/11)")
    print("         = (1/4) * [ 2log(5/4) - log(2) - log(5/11) ]")
    print("         = (1/4) * log[ ((5/4)^2) / (2 * 5/11) ]")
    print("         = (1/4) * log[ (25/16) / (10/11) ]")
    
    print("\nThe final simplified equation is:")
    # Final equation values
    a = 1
    b_eq = 4
    num = 55
    den = 32
    print(f"W* - W = ({a}/{b_eq}) * log({num}/{den})")

    final_result = w_star_val - w_actual_val
    print(f"\nThe numerical value of the difference is approximately: {final_result:.5f}")

solve_race_growth_rate()