import math

def solve_methylation_problem():
    """
    Calculates the initial histone methylation percentage and the resulting
    impact on gene expression based on the problem's parameters.
    """
    print("Part 1: Calculating the initial percentage of trimethylated sites (H3K4me3)")
    print("----------------------------------------------------------------------")

    # Given values for Part 1
    p_final = 11.04  # Percentage of H3K4me3 after 10 hours
    rate = 0.10      # 10% turnover rate per hour
    time = 10        # 10 hours

    # We model the process using continuous exponential decay: P(t) = P0 * e^(-k*t).
    # To find the initial percentage P0, we rearrange the formula to: P0 = P(t) * e^(k*t).
    
    # Calculation
    p_initial = p_final * math.exp(rate * time)
    
    print("The formula for the initial percentage (P0) is: P0 = P(t) * e^(rate * time)")
    print(f"Using the given values:")
    print(f"P0 = {p_final}% * e^({rate} * {time})")
    print(f"P0 = {p_final}% * e^{rate * time}")
    print(f"P0 = {p_final}% * {math.exp(rate * time):.4f}")
    print(f"The initial percentage of trimethylated sites is {p_initial:.2f}%\n")

    # For Part 2, we'll use the clean, rounded value for the initial percentage.
    p_initial_rounded = round(p_initial)
    
    print("Part 2: Calculating the impact on target gene expression")
    print("---------------------------------------------------------")

    # Given values for Part 2
    e_initial = 200      # Initial gene expression in RPKM
    decrease_percentage = 0.10 # 10% relative decrease

    # Calculate the new methylation percentage after a 10% decrease
    p_new = p_initial_rounded * (1 - decrease_percentage)
    
    # Assuming a linear relationship (direct proportionality), E_new = E_initial * (P_new / P_initial)
    e_new = e_initial * (p_new / p_initial_rounded)

    print("Assuming a linear relationship, the new expression level (E_new) is proportional to the new methylation level (P_new).")
    print(f"The initial methylation level (P0) is {p_initial_rounded}%, and initial expression (E0) is {e_initial} RPKM.")
    print(f"After a {decrease_percentage * 100}% relative decrease, the new methylation level is:")
    print(f"P_new = {p_initial_rounded}% * (1 - {decrease_percentage}) = {p_new:.2f}%")
    print(f"\nThe new gene expression is calculated as: E_new = E0 * (P_new / P0)")
    print(f"E_new = {e_initial} * ({p_new:.2f} / {p_initial_rounded})")
    print(f"The new gene expression level is {e_new:.2f} RPKM.")

    # Final answer summary for the <<<>>> tag
    final_answer = f"The initial percentage of trimethylated sites is {p_initial:.2f}%. The new gene expression level would be {e_new:.2f} RPKM."
    return final_answer

# Run the function and capture the final answer string
final_answer_string = solve_methylation_problem()

# The final answer in the required format will be printed by the user's execution environment.
# To conform to the instruction, we'll print it here as well.
# print(f"\n<<<{final_answer_string}>>>")