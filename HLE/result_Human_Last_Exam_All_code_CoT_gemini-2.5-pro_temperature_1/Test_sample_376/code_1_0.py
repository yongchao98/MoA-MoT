def calculate_clonal_expansion():
    """
    Calculates the clonal expansion score based on given genetic data.
    """

    # --- Step 1: Define the given parameters ---

    # Weights for the impact of each gene
    weights = {
        'Oncogene A': 0.5,
        'Tumor suppressor B': -0.7,
        'Oncogene C': 0.4,
        'Tumor suppressor D': -0.6,
        'Oncogene E': 0.3,
        'Tumor suppressor F': -0.5
    }

    # Chromosomal Copy Number Variations (CNVs)
    # Positive for gain, negative for loss
    cnv = {
        'Chromosome 1': 3,  # Gain of 3
        'Chromosome 2': -2, # Loss of 2
        'Chromosome 3': 2   # Gain of 2
    }

    # --- Step 2: Calculate the contribution of each relevant gene ---

    # Oncogene A is on Chromosome 1 (gain of 3 copies)
    copies_gained_A = cnv['Chromosome 1']
    weight_A = weights['Oncogene A']
    score_A = copies_gained_A * weight_A

    # Tumor suppressor B is on Chromosome 2 (loss of 2 copies)
    # Loss of a tumor suppressor promotes expansion, so we flip the sign of the weight.
    copies_lost_B = abs(cnv['Chromosome 2'])
    weight_B = weights['Tumor suppressor B']
    score_B = copies_lost_B * (-weight_B)

    # Oncogene C is on Chromosome 3 (gain of 2 copies)
    copies_gained_C = cnv['Chromosome 3']
    weight_C = weights['Oncogene C']
    score_C = copies_gained_C * weight_C
    
    # Other genes do not contribute to the score as per the rules:
    # - Tumor suppressor D: Gained copies, but weight is for loss (Score = 0)
    # - Oncogene E: Lost copies, but weight is for gain (Score = 0)
    # - Tumor suppressor F: Gained copies, but weight is for loss (Score = 0)

    # --- Step 3: Sum the scores and print the results ---

    total_score = score_A + score_B + score_C

    print("Clonal Expansion Score Calculation:")
    print("The score is the sum of impacts from Oncogene A (gain), Tumor Suppressor B (loss), and Oncogene C (gain).")
    print("Other genes have no impact as their changes (gain/loss) do not match their weight definitions.\n")

    print("Final Equation:")
    # Print each number in the equation
    print(f"Score = (copies_A * weight_A) + (copies_B_lost * -weight_B) + (copies_C * weight_C)")
    print(f"Score = ({copies_gained_A} * {weight_A}) + ({copies_lost_B} * {-weight_B}) + ({copies_gained_C} * {weight_C})")
    print(f"Score = {score_A} + {score_B} + {score_C}")
    print(f"\nTotal Clonal Expansion Score = {total_score}")

calculate_clonal_expansion()
<<<3.7>>>