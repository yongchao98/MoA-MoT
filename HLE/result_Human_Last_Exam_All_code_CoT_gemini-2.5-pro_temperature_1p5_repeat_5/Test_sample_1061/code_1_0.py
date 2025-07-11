import textwrap

def analyze_tree_ring_question():
    """
    Analyzes a multiple-choice question about tree ring isotopes
    and determines the most likely answer.
    """
    question = "Which of the factors below was the predominant factor influencing the declining 13C ratio observed in the tree ring data over this period (1886-1990 AD)?"
    
    options = {
        'A': 'An increase in tree ring thickness as the tree matures',
        'B': 'Periods of drought',
        'C': 'Increased photosynthetic reserves of starch fueling tree growth',
        'D': 'Thinning earlywood tree ring proportion',
        'E': 'Changes in the SE Asia monsoon'
    }

    print("Analyzing the question: '{}'".format(question))
    print("-" * 30)
    print("Core Observation: A declining 13C ratio in Chinese pine trees from 1886-1990.\n")
    
    print("Key Scientific Principles:")
    print("1. The 'Suess Effect': Burning fossil fuels since the Industrial Revolution has lowered the 13C ratio in the entire atmosphere. Trees record this global trend.")
    print("2. Climate Influence: Water availability affects how trees take up carbon. More water (e.g., from a strong monsoon) leads to a lower 13C ratio. Drought leads to a higher 13C ratio.\n")
    
    print("Evaluating the Options:")
    print("-" * 30)

    # Analysis of each option
    reasoning = {
        'A': "This is a biological effect specific to an individual tree's age. It is not a predominant factor for a widespread, century-long regional trend.",
        'B': "Drought causes stomata to close, which INCREASES the 13C ratio. This is the opposite of the observed decline.",
        'C': "This is an internal physiological factor and is unlikely to be the primary driver of a long-term, large-scale trend.",
        'D': "Changes in wood type proportions can affect the 13C ratio, but this is a secondary effect and not the predominant driver of a major regional trend.",
        'E': "The SE Asia monsoon is the dominant climatic driver of precipitation in this region. Changes in monsoon strength directly impact water availability, which in turn is a primary control on the 13C ratio in tree rings (stronger monsoon = lower 13C). While the 'Suess Effect' is the ultimate cause of the baseline decline, it is not an option. Of the choices given, the monsoon is the most significant environmental factor influencing the 13C record."
    }

    eliminated_options = []
    best_option = None
    
    # Option B is factually opposite
    print(f"Option B analysis: '{options['B']}' is incorrect. {reasoning['B']}\n")
    eliminated_options.append('B')
    
    # Options A, C, D are less significant
    for opt in ['A', 'C', 'D']:
        print(f"Option {opt} analysis: '{options[opt]}' is a minor factor. {reasoning[opt]}\n")
        eliminated_options.append(opt)

    # Option E is the best fit
    best_option = 'E'
    print(f"Option E analysis: '{options['E']}' is the most plausible answer. {reasoning['E']}\n")

    print("-" * 30)
    print("Conclusion:")
    print("Options A, B, C, and D are either incorrect or represent minor factors.")
    print("Option E represents the major regional climate system that directly influences the conditions controlling the 13C ratio in trees.")
    print("Therefore, it is the best answer among the choices provided.")
    print("-" * 30)
    print(f"Final Answer Choice: {best_option}")

if __name__ == '__main__':
    analyze_tree_ring_question()