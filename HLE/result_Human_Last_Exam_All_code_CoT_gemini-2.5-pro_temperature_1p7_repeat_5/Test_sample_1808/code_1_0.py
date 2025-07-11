def analyze_hybrid_zone_genetics():
    """
    Analyzes the effects of gene flow in a hybrid zone on various genetic metrics.
    """
    print("Analyzing the potential outcomes of gene flow across a hybrid zone:")
    print("="*60)

    # Dictionary of explanations
    explanations = {
        "A. High Fst": "Fst measures population differentiation. While gene flow reduces Fst, strong selection against hybrids can maintain high differentiation. So, this CAN occur.",
        "B. High Dxy": "Dxy measures absolute divergence, reflecting historical separation. Gene flow does not erase this history. So, this CAN occur.",
        "C. High Fis": "Fis measures heterozygote deficit. The Wahlund effect in a mixed sample from a hybrid zone leads to a high Fis. So, this CAN occur.",
        "D. High u (mutation rate)": "The mutation rate (u) is a fundamental biochemical property. Gene flow moves existing alleles; it does not alter the rate of new mutations. So, this CANNOT occur as a result of gene flow.",
        "E. High Pi": "Pi measures genetic diversity within a population. Gene flow introduces new alleles, increasing Pi. So, this CAN occur."
    }

    correct_answer = None
    for option, reason in explanations.items():
        print(f"Option {option}:\n   - {reason}\n")
        if "CANNOT occur" in reason:
            correct_answer = option

    print("="*60)
    print(f"Conclusion: The event that cannot occur as a result of gene flow is a high mutation rate.")
    print(f"The correct option is: {correct_answer.split('.')[0]}")


if __name__ == "__main__":
    analyze_hybrid_zone_genetics()