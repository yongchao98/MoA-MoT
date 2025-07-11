import textwrap

def find_common_element():
    """
    Analyzes the cinematic elements shared between Fritz Lang and William Friedkin.
    """
    director_1 = "Fritz Lang"
    director_2 = "William Friedkin"
    
    analysis = {
        'A': "Incorrect. Neither director is known for using Aboriginal masks.",
        'B': "Incorrect. Magic wands are not a feature of either director's style.",
        'C': "Incorrect. While Lang's 'Metropolis' features a famous cinematic robot, Friedkin's films do not contain cyborgs.",
        'D': "Correct. Friedkin directed 'Bug' (2006), about a delusion of insects. Lang's 'The Testament of Dr. Mabuse' (1933) features a famous scene where a character's madness is visualized as a tormenting swarm of creatures, often interpreted as bugs. This represents a shared visual motif for psychological horror.",
        'E': "Incorrect, as a common element was found."
    }
    
    correct_option = 'D'
    explanation = analysis[correct_option]
    
    print(f"Analyzing the works of {director_1} and {director_2}...\n")
    print(textwrap.fill(f"The correct answer is {correct_option}. {explanation}", width=80))

find_common_element()