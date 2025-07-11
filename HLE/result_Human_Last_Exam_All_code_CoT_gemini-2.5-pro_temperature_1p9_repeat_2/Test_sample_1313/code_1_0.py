def find_restored_scene():
    """
    Analyzes the provided information to determine the content of the
    newly restored scene in Fritz Lang's "Kriemhild's Revenge".
    """
    # The key information is in the Le Monde article [2], which describes the restored version.
    # The article contrasts the previously known ending with the new one.
    
    # Old Ending Description (Paraphrased from article):
    # King Etzel is seen alone amidst the carnage, lamenting the death of his wife.
    
    # New Ending Description (Paraphrased from article):
    # The restored ending shows Etzel picking up Kriemhild's crown and giving it to Hildebrand.
    # Crucially, it then shows him carrying the body of his dead son, Ortlieb, in his arms.
    # He then renounces his throne.
    
    print("Analysis of the provided sources reveals the following:")
    print("The Le Monde article about the 2010 restoration explicitly describes the newly added sequence at the end of the film.")
    print("According to the article, after Kriemhild's death, the restored footage shows King Etzel taking several actions.")
    
    # Matching the new ending elements to the answer choices:
    # A. A shot of Hildebrand bringing Kriemhild's crown to Etzel. -> Incorrect. Etzel gives the crown to Hildebrand.
    # B. A shot of Etzel watching the sunset, mired in sorrow. -> Incorrect. No sunset is mentioned.
    # C. A shot of Hildebrand striking Kriemhild down with his spear. -> Incorrect. This event precedes the final scene and is part of the original narrative.
    # D. A shot of Etzel lifts his infant son amidst the carnage. -> Correct. The article states he carries the body of his son, Ortlieb, in his arms.
    # E. A shot of Etzel calling for help, lamenting. -> Incorrect. This was part of the previously known ending.
    # F. A shot of Kriemhild placing Siegfried's sword on her lap in the shadowy hall. -> Incorrect. Not mentioned in the description of the ending.

    correct_choice = "D"
    
    print(f"One of the key elements in this new sequence is Etzel carrying the body of his son.")
    print(f"This directly corresponds to choice '{correct_choice}'.")

find_restored_scene()