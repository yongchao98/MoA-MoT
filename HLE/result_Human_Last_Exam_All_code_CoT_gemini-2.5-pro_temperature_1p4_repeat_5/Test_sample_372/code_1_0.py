import math

def solve_epigenetics_problem():
    """
    Solves a two-part problem related to histone methylation and gene expression.
    """
    
    # Part 1: Calculate the initial percentage of H3K4me3 sites.
    # Given values for Part 1
    p_t_percent = 11.04  # Percentage of H3K4me3 sites after 10 hours
    decay_rate = 0.10      # 10% decay per hour
    time = 10              # 10 hours
    
    # Convert percentage to a proportion for calculation
    p_t_proportion = p_t_percent / 100
    
    # Calculate the initial proportion using the exponential decay formula: P(t) = P(0) * (1 - r)^t
    # Rearranged formula: P(0) = P(t) / (1 - r)^t
    initial_proportion = p_t_proportion / math.pow(1 - decay_rate, time)
    initial_percentage = initial_proportion * 100
    
    print("--- Part 1: Initial Percentage of H3K4me3 ---")
    print("The initial percentage of trimethylated sites is calculated using the exponential decay formula.")
    print("Calculation:")
    print(f"Initial Percentage = Final Percentage / (1 - Decay Rate)^Time")
    print(f"Initial Percentage = {p_t_percent}% / (1 - {decay_rate})^{time} = {initial_percentage:.2f}%")
    print("-" * 50)
    
    # Part 2: Determine the impact on gene expression.
    # Given values for Part 2
    initial_expression = 200.0  # RPKM
    proportion_decrease = 0.10    # 10% decrease
    
    # A linear relationship implies that a 10% relative decrease in the methylation
    # proportion will lead to a 10% relative decrease in gene expression.
    new_expression = initial_expression * (1 - proportion_decrease)
    
    print("--- Part 2: Impact on Gene Expression ---")
    print("Due to the linear relationship, a 10% decrease in the proportion of H3K4me3 sites")
    print("results in a corresponding 10% decrease in gene expression.")
    print("Calculation:")
    print(f"New Expression = Initial Expression * (1 - Proportional Decrease)")
    print(f"New Expression = {initial_expression} RPKM * (1 - {proportion_decrease}) = {new_expression:.1f} RPKM")
    
    # Final answer based on Part 2
    final_answer = new_expression
    print(f"\n<<<>>>\n") # An empty line is added here for better formatting before the final answer.
    print(f"<<<{final_answer}>>>")

solve_epigenetics_problem()