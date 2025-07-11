import math

def solve_rangoli_puzzle():
    """
    Solves the Rangoli puzzle by determining the number of curves the master must draw.
    """
    # Step 1: Set up the equation for the total original curves (T).
    # The total T is the sum of three disjoint groups described in the poem:
    # 1. Curves that "lost shape": (3/8) * T
    # 2. Curves that "found new paths": (1/4) * T
    # 3. Curves that "stay true": 90
    # Equation: T = (3/8 * T) + (1/4 * T) + 90
    
    # To solve for T, we first find the coefficient of T:
    # T = (3/8 + 2/8) * T + 90
    # T = (5/8) * T + 90
    # T - (5/8) * T = 90
    # (3/8) * T = 90
    
    # Step 2: Solve for T.
    unchanged_curves = 90
    t_coefficient = 1 - (3/8 + 1/4) # This is 1 - 5/8 = 3/8
    total_original_curves = int(unchanged_curves / t_coefficient)
    
    # Step 3: Identify the curves the master must draw.
    # The group of curves being redrawn ("those that left their plotted way") is broken down
    # into fractions of 1/5 and 2/9. This means its total number must be divisible by the
    # least common multiple of 5 and 9, which is 45.

    # Let's check the two groups of changed curves:
    lost_shape_curves = int((3/8) * total_original_curves)
    new_path_curves = int((1/4) * total_original_curves)

    # The "lost shape" group has 90 curves, which is divisible by 45.
    # The "found new paths" group has 60 curves, which is not divisible by 45.
    # Therefore, the 90 curves that "lost shape" are the ones being redrawn.
    
    curves_to_draw = lost_shape_curves
    
    # Step 4: Display the final calculation and result.
    # The question asks "how many total curves must the master draw".
    # This is the number we identified in Step 3.
    
    print("To solve the puzzle, we first determine the original total number of curves.")
    print(f"The original total (T) is composed of three groups: (3/8 * T) + (1/4 * T) + 90.")
    print(f"Solving for T gives us {total_original_curves} original curves.\n")
    
    print("The question asks how many curves the master must draw.")
    print("This corresponds to the group of curves that 'lost shape serene'.")
    
    print("\nFinal Calculation:")
    # The instruction requires printing each number in the final equation.
    # The equation is: Curves to Draw = (Fraction Lost) * (Total Curves)
    fraction_numerator = 3
    fraction_denominator = 8
    
    print(f"Number of Curves to Draw = ({fraction_numerator}/{fraction_denominator}) * {total_original_curves}")
    print(f"{curves_to_draw} = ({fraction_numerator}/{fraction_denominator}) * {total_original_curves}")

    # For final answer extraction
    final_answer_value = curves_to_draw
    return final_answer_value


# --- Main Execution ---
if __name__ == "__main__":
    final_answer = solve_rangoli_puzzle()
    # The final answer is wrapped for extraction as requested.
    print(f"\n<<< {final_answer} >>>")