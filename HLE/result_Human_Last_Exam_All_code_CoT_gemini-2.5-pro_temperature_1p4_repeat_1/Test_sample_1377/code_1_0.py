import math

def solve_race_problem():
    """
    Calculates the difference between optimal and actual investment growth rates
    based on true vs. mistaken probabilities in a three-competitor race.
    """
    # Step 1: Define problem parameters
    # True probabilities for competitors (A, B, C)
    p_true = [1/2, 1/4, 1/4]
    # Mistakenly believed probabilities for (A, B, C)
    p_believed = [1/4, 1/2, 1/4]
    # Net odds (payout is b:1, meaning you get b+1 back for a 1 unit bet)
    b = [3, 2, 2]
    
    print("This program calculates W* - W, the difference between the optimal and actual growth rates.")
    print("The final equation is of the form: (a/b)*ln(c) - (d/e)*ln(f)\n")

    # --- Part 1: Calculate the optimal strategy and growth rate (W*) ---
    
    # Under true probabilities, an edge exists if p_true * (b+1) > 1
    # For A: (1/2) * (3+1) = 2.0 > 1. Positive edge.
    # For B: (1/4) * (2+1) = 0.75 <= 1. No edge.
    # For C: (1/4) * (2+1) = 0.75 <= 1. No edge.
    # Therefore, the optimal strategy is to only bet on A.

    # The Kelly fraction for a single bet is f = p - (1-p)/b
    f_A_star = p_true[0] - (1 - p_true[0]) / b[0]
    
    # Optimal growth rate W* = p_A*ln(1+b_A*f_A*) + (p_B+p_C)*ln(1-f_A*)
    w_star = (p_true[0] * math.log(1 + b[0] * f_A_star) +
              (p_true[1] + p_true[2]) * math.log(1 - f_A_star))

    # --- Part 2: Calculate the mistaken strategy and actual growth rate (W) ---

    # Under believed probabilities, an edge exists if p_believed * (b+1) > 1
    # For A: (1/4) * (3+1) = 1.0 <= 1. No edge.
    # For B: (1/2) * (2+1) = 1.5 > 1. Positive edge.
    # For C: (1/4) * (2+1) = 0.75 <= 1. No edge.
    # The mistaken strategy is to only bet on B.
    
    f_B_mistaken = p_believed[1] - (1 - p_believed[1]) / b[1]
    
    # The actual growth rate W is calculated using mistaken fractions but with true probabilities.
    # Wealth factor if A or C wins (bet on B is lost): 1 - f_B_mistaken
    # Wealth factor if B wins: 1 + b_B * f_B_mistaken
    w = (p_true[0] * math.log(1 - f_B_mistaken) +
         p_true[1] * math.log(1 + b[1] * f_B_mistaken) +
         p_true[2] * math.log(1 - f_B_mistaken))

    # --- Part 3: Compute and display the difference ---
    
    difference = w_star - w
    
    # The symbolic form of the difference is (11/4)*ln(2) - (3/2)*ln(3)
    c1, n1, d1 = 11, 4, 2
    c2, n2, d2 = 3, 2, 3

    print("Step-by-step derivation:")
    print(f"Optimal growth rate W* = ln(2) - (1/2)*ln(3) ≈ {w_star:.5f}")
    print(f"Actual growth rate W = ln(3) - (7/4)*ln(2) ≈ {w:.5f}")
    print("\nThe difference W* - W is derived as follows:")
    print(f"(ln(2) - (1/2)*ln(3)) - (ln(3) - (7/4)*ln(2))")
    print(f"= (1 + 7/4)*ln(2) - (1/2 + 1)*ln(3)")
    print(f"= ({c1}/{n1})*ln({d1}) - ({c2}/{n2})*ln({d2})")
    
    print("\nFinal equation and result:")
    print(f"({c1}/{n1}) * ln({d1}) - ({c2}/{n2}) * ln({d2}) = {difference:.5f}")
    
    # Final answer in requested format
    return f"<<<{difference:.5f}>>>"

# Execute the function and print the final answer
result_string = solve_race_problem()
print(result_string)
