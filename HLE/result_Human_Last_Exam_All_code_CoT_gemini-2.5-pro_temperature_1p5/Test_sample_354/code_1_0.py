def solve_sociology_question():
    """
    Analyzes Alan Bryman's concept of Disneyization to answer the multiple-choice question.
    """

    # Step 1: Define Alan Bryman's four key dimensions of Disneyization.
    bryman_dimensions = {
        "Theming": "The application of a narrative to a place or object.",
        "Hybrid Consumption": "The blurring of boundaries between different forms of consumption (e.g., a shop is also an attraction).",
        "Merchandising": "The promotion and sale of goods bearing logos and imagery.",
        "Performative Labor": "Making employees part of the performance through scripts, costumes, and controlled emotional displays."
    }

    # Step 2: Analyze the provided answer choices.
    answer_choices = {
        "A": ["hybrid consumption", "merchandising"],
        "B": ["performative labor", "sanitization"],
        "C": ["trivialization", "theming"],
        "D": ["sanitization", "trivialization"],
        "E": ["disneyfication", "disneyization"],
        "F": ["mcdonaldization", "disneyization"],
        "G": ["theming", "performative labor"]
    }

    # Step 3: Identify the best fit.
    # We are looking for an option where both items are in Bryman's list of core dimensions.
    # Choice A, "hybrid consumption and merchandising," are both correct dimensions.
    # Choice G, "theming and performative labor," are also both correct dimensions.
    # However, theming (creating the environment) and performative labor (populating it with actors)
    # are often considered the most central and defining characteristics of the process.
    
    # Step 4: Select the chosen characteristics and final answer.
    chosen_characteristic_1 = "theming"
    chosen_characteristic_2 = "performative labor"
    final_answer = "G"

    # Step 5: Print the explanation and the result.
    print(f"In his book 'The Disneyization of Society' (2004), Alan Bryman identifies four main characteristics of Disneyization.")
    print(f"Two of these core characteristics are {chosen_characteristic_1} and {chosen_characteristic_2}.")
    print(f"Theming refers to imposing a narrative on a space, while performative labor involves scripting employee behavior to be part of the show.")
    print("Based on this, the correct choice is G.")

solve_sociology_question()
print("\n<<<G>>>")