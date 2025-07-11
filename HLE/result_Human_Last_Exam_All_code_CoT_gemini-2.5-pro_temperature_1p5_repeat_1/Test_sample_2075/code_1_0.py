import pandas as pd
import io

def find_purest_emotion_region():
    """
    This function identifies the dmPFC sub-region most purely activated by emotion
    by using a conceptual model based on established neuroscience findings.
    """
    # Step 1: Create a knowledge base representing brain region functions.
    # The 'emotion_purity_score' is a conceptual value from 0.0 to 1.0.
    # It represents how selectively the region is involved in emotion processing.
    # A high score (e.g., 0.9) means the region's primary, most specific role is emotion.
    # A moderate score (e.g., 0.5) means emotion is a key function, but heavily mixed with others (like cognition).
    # A low score (e.g., 0.1) means the region is not primarily involved in emotion.
    data = """
region,primary_function,emotion_purity_score
caudal-right,action_monitoring,0.1
rostroventral,emotion_valuation,0.9
rostrodorsal,social_cognition,0.5
caudal-left,action_monitoring,0.1
occipital,vision,0.0
"""
    # Using io.StringIO to simulate reading from a file
    knowledge_base = pd.read_csv(io.StringIO(data))

    # Step 2: Find the region with the maximum 'emotion_purity_score'.
    most_pure_region_row = knowledge_base.loc[knowledge_base['emotion_purity_score'].idxmax()]
    best_region_name = most_pure_region_row['region']

    # Step 3: Print the reasoning and the result.
    print("Analysis based on functional specialization of dmPFC sub-regions:")
    print("-" * 60)
    for index, row in knowledge_base.iterrows():
        print(f"Region: {row['region']:<15} | Primary Function: {row['primary_function']:<20} | Emotion Purity Score: {row['emotion_purity_score']}")
    print("-" * 60)
    print(f"\nThe question asks for the sub-region most *purely* activated by emotion processing.")
    print(f"Based on our knowledge model, the region with the highest emotion purity score is '{best_region_name}'.")
    print(f"\nConclusion: The {best_region_name} dmPFC is most consistently and selectively associated with core emotion and valuation processes, making it the area most purely activated by such stimuli.")

find_purest_emotion_region()
<<<B>>>