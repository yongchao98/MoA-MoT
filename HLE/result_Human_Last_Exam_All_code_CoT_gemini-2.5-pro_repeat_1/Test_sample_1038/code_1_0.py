import pandas as pd

def evaluate_design_strategies():
    """
    Analyzes and scores different circuit design strategies based on a set of criteria.
    """
    # Define the design options and the evaluation criteria
    options = {
        'A': "Use minimum-length transistors with large width and a small gate capacitor.",
        'B': "Split the gate capacitor into refreshed segments.",
        'C': "Use an on-chip body-bias generator to increase threshold voltage.",
        'D': "Replace bootstrapping with high-impedance current mirrors.",
        'E': "Use a split-gate transistor with mixed static/bootstrapped control."
    }
    
    criteria = ['Headroom', 'Offset Tolerance', 'Reset Time', 'Leakage Control']
    
    # Scoring based on the analysis: -2 (Very Poor) to +2 (Excellent)
    scores = {
        'A': [-1, -1, 2, -2],
        'B': [0, -1, 0, 1],
        'C': [2, 1, 0, 1],
        'D': [-2, -2, 1, 2],
        'E': [-1, -1, 1, -2]
    }
    
    # Weights for each criterion. Core operational metrics are weighted higher.
    weights = {
        'Headroom': 2,
        'Offset Tolerance': 2,
        'Reset Time': 1,
        'Leakage Control': 2
    }
    
    print("Evaluating design strategies for the bootstrapped pseudo-resistor...\n")
    
    results = {}
    best_option = None
    max_score = -float('inf')

    # Calculate and display scores
    for option_key, option_desc in options.items():
        print(f"--- Option {option_key}: {option_desc} ---")
        
        total_score = 0
        equation_parts = []
        
        for i, criterion in enumerate(criteria):
            score = scores[option_key][i]
            weight = weights[criterion]
            total_score += score * weight
            
            # Formatting the equation part to show numbers
            part_str = f"({criterion}: {score} * W:{weight})"
            equation_parts.append(part_str)

        equation_str = " + ".join(equation_parts)
        print(f"Calculation: {equation_str} = {total_score}")
        results[option_key] = total_score
        
        if total_score > max_score:
            max_score = total_score
            best_option = option_key
        print("-" * (len(option_desc) + 18) + "\n")


    print("--- Conclusion ---")
    print("The goal is to find the strategy that most effectively balances all requirements.")
    print("Option C (Body-Bias Generator) receives the highest score. Here's why:")
    print("- It directly improves headroom and leakage control by increasing the threshold voltage (Vt).")
    print("- A higher Vt makes the device more tolerant to the Vgs variations caused by sensor offsets.")
    print("- While it adds complexity, this approach provides the best balance for the critical operational parameters at a low 1.2V supply.")

    print(f"\nThe most effective strategy is Option {best_option}.")


if __name__ == "__main__":
    evaluate_design_strategies()
