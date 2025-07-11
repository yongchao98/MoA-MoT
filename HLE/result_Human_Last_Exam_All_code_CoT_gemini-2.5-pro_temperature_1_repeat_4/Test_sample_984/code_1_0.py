import math

def analyze_sorting_efficiency():
    """
    Analyzes and explains the low efficiency of a cell sorting experiment
    based on the principle of rare biological interactions vs. common random associations.
    """

    # --- Step 1: Define the experimental results from the problem ---
    sorted_wells = 95
    successful_wells = 4

    # --- Step 2: Calculate the observed experimental efficiency ---
    experimental_efficiency = successful_wells / sorted_wells

    print("--- Analysis of Cell Sorting Experiment ---")
    print(f"\n1. Calculating the observed experimental efficiency:")
    print(f"   Successfully recovered paired chains from {successful_wells} wells.")
    print(f"   Total wells sorted with doublets: {sorted_wells}.")
    print(f"   Equation: Efficiency = Successful Wells / Total Sorted Wells")
    print(f"   Calculation: Efficiency = {successful_wells} / {sorted_wells} = {experimental_efficiency:.4f}")
    print(f"   This means the success rate was {experimental_efficiency:.2%}.")

    # --- Step 3: Propose a hypothesis based on Choice B ---
    print("\n2. Hypothesis: True antigen-specific interactions are rare events, and many observed doublets are random associations.")
    print("   The sorter cannot distinguish a true interaction from a random collision.")
    print("   Therefore, the final efficiency is determined by the proportion of true doublets among all sorted doublets.")

    # --- Step 4: Calculate the ratio of random vs. true doublets from the data ---
    # Efficiency = True / (True + Random)
    # If we consider the number of True events as 1, then the number of Random events is R.
    # Efficiency = 1 / (1 + R) --> R = (1 / Efficiency) - 1
    if experimental_efficiency > 0:
        random_to_true_ratio = (1 / experimental_efficiency) - 1
    else:
        random_to_true_ratio = float('inf')

    print("\n3. Estimating the prevalence of random vs. true doublets:")
    print(f"   Let R be the number of random doublets for every 1 true doublet.")
    print(f"   Equation: Efficiency = 1 / (1 + R)")
    print(f"   Solving for R: R = (1 / Efficiency) - 1")
    print(f"   Calculation: R = (1 / {experimental_efficiency:.4f}) - 1 = {random_to_true_ratio:.2f}")
    print(f"\n   This result suggests that for every 1 true, specific doublet sorted,")
    print(f"   the machine also sorted approximately {math.ceil(random_to_true_ratio)} random, non-specific doublets.")

    # --- Step 5: Final Conclusion ---
    print("\n4. Conclusion:")
    print("   The low efficiency is not primarily due to technical failure but is a direct consequence of the biology.")
    print("   The vast majority of events identified as 'doublets' by the sorter were likely random cell pairings.")
    print("   The 4% success rate accurately reflects the rarity of the true antigen-specific B cell:tumor cell interactions within the total population of sorted events.")
    print("   This strongly supports answer choice B.")

if __name__ == '__main__':
    analyze_sorting_efficiency()