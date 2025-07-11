import pandas as pd

def find_contrasting_feature():
    """
    This function models the effects of various linguistic features on word acquisition
    for children versus language models (LLMs) to identify the one that fits the problem description.
    """
    # Data representing the known effects of word features on acquisition.
    # +1: Strong Positive Effect (easier to learn)
    # -1: Strong Negative Effect (harder to learn)
    #  0: Neutral, mixed, or effect that is not oppositional (e.g., positive for both)
    data = {
        'Feature': [
            'word concreteness rating',
            'number of characters',
            'mean length of sequences containing a given word',
            'word complexity rating via Flesch Reading Ease',
            'word imageability rating',
            'lexical class',
            'word log frequency'
        ],
        'Answer Choice': ['A', 'B', 'C', 'D', 'E', 'F', 'G'],
        # Effect on Children: Concrete/imageable/frequent words are easier (+). Long/complex words are harder (-).
        'Effect on Child Acquisition': [1, -1, -1, -1, 1, 0, 1],
        # Effect on LLMs: Frequent/context-rich words are easier (+). Context-variable (concrete) words are harder (-).
        'Effect on LLM Acquisition': [-1, -1, 1, -1, -1, 0, 1]
    }

    df = pd.DataFrame(data)

    # The problem asks for a feature with:
    # 1. A strong significant POSITIVE effect on children.
    # 2. A strong significant NEGATIVE effect on language models.
    
    target_condition = (df['Effect on Child Acquisition'] > 0) & (df['Effect on LLM Acquisition'] < 0)
    
    result = df[target_condition]

    print("Analyzing features based on findings from psycholinguistics and machine learning...")
    print("="*80)
    print("Searching for a feature with a POSITIVE effect on children AND a NEGATIVE effect on LLMs.")
    print("="*80)

    if not result.empty:
        print("Found matching feature(s):")
        for index, row in result.iterrows():
            print(f"\nAnswer Choice: {row['Answer Choice']}")
            print(f"Feature Name: {row['Feature']}")
            print("  - Explanation:")
            print("    - This feature has a POSITIVE effect on child acquisition because children learn tangible,")
            print("      easily visualized concepts more quickly.")
            print("    - This feature has a NEGATIVE effect on LLM acquisition because such words appear in highly")
            print("      diverse contexts, making them harder for a model to learn by predicting from surrounding text.")

    else:
        print("No feature matching the criteria was found based on the modeled effects.")

find_contrasting_feature()