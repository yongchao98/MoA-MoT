import math

def solve_methylation_expression():
    """
    Solves a two-part problem regarding histone methylation and gene expression.
    Part 1: Calculates the initial percentage of H3K4me3 sites based on a decay model.
    Part 2: Calculates the new gene expression level based on a decrease in H3K4me3 sites.
    """

    # --- Part 1: Calculate initial methylation percentage ---

    print("--- Part 1: Initial Trimethylation Percentage ---")
    
    # Given parameters
    final_proportion_percent = 11.04
    decay_rate_per_hour = 0.10  # 10%
    time_hours = 10

    # The decay model is P_final = P_initial * (1 - rate)^time
    # We rearrange to find P_initial = P_final / (1 - rate)^time
    
    # Calculate the decay factor (1 - rate)^time
    decay_factor = (1 - decay_rate_per_hour)**time_hours
    
    # Calculate the initial proportion in percent
    initial_proportion_percent = final_proportion_percent / decay_factor

    print("The final percentage of H3K4me3 sites is given as 11.04% after 10 hours with a decay rate of 10% per hour.")
    print("The calculation for the initial percentage is:")
    # Print the equation with all the numbers
    print(f"Initial Percentage = {final_proportion_percent} / (1 - {decay_rate_per_hour}) ** {time_hours}")
    print(f"Initial Percentage = {final_proportion_percent} / {decay_factor:.4f}")
    print(f"The initial percentage of trimethylated sites was: {initial_proportion_percent:.2f}%\n")


    # --- Part 2: Impact on gene expression ---
    
    print("--- Part 2: Impact on Gene Expression ---")

    # Given parameters
    initial_expression_rpkm = 200
    proportion_decrease_percent = 10 # 10%

    # Calculate the new expression level
    # New Expression = Initial Expression * (1 - decrease_fraction)
    decrease_fraction = proportion_decrease_percent / 100.0
    new_expression_rpkm = initial_expression_rpkm * (1 - decrease_fraction)
    
    print("Assuming a linear relationship, a 10% decrease in H3K4me3 sites causes a 10% decrease in gene expression.")
    print("The initial expression is 200 RPKM.")
    print("The calculation for the new expression level is:")
    # Print the equation with all the numbers
    print(f"New Expression = {initial_expression_rpkm} * (1 - {decrease_fraction})")
    print(f"The new gene expression level is: {new_expression_rpkm:.0f} RPKM.")


solve_methylation_expression()
# The final answer contains the result from both parts of the question.
final_answer_text = "The initial percentage of trimethylated sites was 31.66%, and the new gene expression level is 180 RPKM."
print(f"\n<<<{final_answer_text}>>>")