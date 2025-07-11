import sys

def solve_color_question():
    """
    This script determines the color of a concentrated Pseudomonas aeruginosa sample
    by explaining the biology of its pigment production.
    """
    
    # Step 1: Identify the key organism and its properties.
    organism = "Pseudomonas aeruginosa"
    print(f"Organism in question: {organism}")
    print("This bacterium is known for producing water-soluble pigments.")
    print("-" * 20)

    # Step 2: Define the primary pigments and their colors.
    # We will use numbers symbolically for the "equation" part of the request.
    pigment1_name = "Pyocyanin"
    pigment1_color = "Blue"
    pigment1_code = 1
    
    pigment2_name = "Pyoverdine"
    pigment2_color = "Yellow-Green"
    pigment2_code = 2
    
    print(f"It produces two main pigments:")
    print(f" - {pigment1_name}, which has a {pigment1_color} color.")
    print(f" - {pigment2_name}, which has a {pigment2_color} color.")
    print("-" * 20)
    
    # Step 3: Explain how the pigments combine to produce the final color.
    # The washing steps remove the medium, but the dense cell pellet
    # retains significant color from the combination of pigments.
    print("When these two pigments are produced together by a dense culture, they mix.")
    print("The combination of blue and yellow-green results in a blue-green appearance.")
    
    # Step 4: Formulate and display the final "equation" as requested.
    result_code = pigment1_code + pigment2_code
    print("\nSymbolic color equation:")
    print(f"{pigment1_code} ({pigment1_color}) + {pigment2_code} ({pigment2_color}) = {result_code} (Blue-Green)")
    print("-" * 20)

    # Step 5: State the final answer based on the biological reasoning.
    final_answer = "C. Blue-green"
    print(f"Therefore, the color of the dense, concentrated sample is Blue-green.")
    print(f"The correct answer choice is: {final_answer}")

solve_color_question()
<<<C>>>