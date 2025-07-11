import math

def solve_epigenetics_problem():
    """
    Solves a two-part problem related to histone methylation and gene expression.
    """
    
    # --- Part 1: Calculate initial percentage of trimethylated sites ---
    print("--- Part 1: Calculating Initial Percentage of Trimethylated Sites ---")

    p_final_percent = 11.04
    time_hours = 10
    decay_rate_per_hour = 0.10

    # Convert final percentage to a proportion for calculation
    p_final_proportion = p_final_percent / 100.0

    # The decay model is P(t) = P(0) * (1 - rate)^t
    # We solve for P(0): P(0) = P(t) / (1 - rate)^t
    
    decay_factor = (1 - decay_rate_per_hour) ** time_hours
    p_initial_proportion = p_final_proportion / decay_factor
    p_initial_percent = p_initial_proportion * 100.0

    print("The final proportion of H3K4me3 after {} hours is {}%.".format(time_hours, p_final_percent))
    print("The rate of conversion to H3K4me2 is {}% per hour.".format(decay_rate_per_hour * 100))
    print("\nThe equation to find the initial proportion P(0) is: P(0) = P(t) / (1 - rate)^t")
    print("Plugging in the values:")
    print("P(0) = {} / (1 - {})**{}".format(p_final_proportion, decay_rate_per_hour, time_hours))
    print("P(0) = {} / {}".format(p_final_proportion, decay_factor))
    print("\nThe calculated initial proportion is {:.4f}, which is {:.2f}%.".format(p_initial_proportion, p_initial_percent))

    print("\n" + "="*70 + "\n")

    # --- Part 2: Determine the impact on gene expression ---
    print("--- Part 2: Determining the Impact on Gene Expression ---")

    initial_expression_rpkm = 200
    proportion_decrease = 0.10  # 10% decrease

    # Due to the linear relationship, a 10% decrease in methylation proportion
    # results in a 10% decrease in gene expression.
    final_expression_rpkm = initial_expression_rpkm * (1 - proportion_decrease)

    print("The initial gene expression level is {} RPKM.".format(initial_expression_rpkm))
    print("The proportion of H3K4me3 sites decreases by {}%.".format(proportion_decrease * 100))
    print("\nGiven the linear relationship, the new expression level is calculated as:")
    print("New Expression = Initial Expression * (1 - Decrease)")
    print("New Expression = {} * (1 - {})".format(initial_expression_rpkm, proportion_decrease))
    print("New Expression = {} * {}".format(initial_expression_rpkm, (1 - proportion_decrease)))
    
    print("\nThe resulting gene expression level is {:.0f} RPKM.".format(final_expression_rpkm))
    
    # Return the final answer for the 'impact' question
    return final_expression_rpkm

# Execute the function and capture the final answer for the 'answer' block
final_answer = solve_epigenetics_problem()