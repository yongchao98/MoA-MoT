import sys

def analyze_covalency():
    """
    Analyzes the relative covalency of two cerium compounds based on orbital overlap.
    """
    # The problem states that CeF6(2-) has greater 4f-2p overlap than CeCl6(2-) has 4f-3p overlap.
    # We can represent this relationship with numerical values. The actual values are not important,
    # only their relative magnitude.
    orbital_overlap_CeF6 = 1.0  # Represents greater overlap
    orbital_overlap_CeCl6 = 0.5  # Represents weaker overlap

    # The strength of covalency is directly proportional to the extent of orbital overlap.
    # We can use an if-statement to compare them and draw a conclusion.

    print("Comparing the covalency based on orbital overlap:")
    print(f"Relative overlap in CeF6(2-): {orbital_overlap_CeF6}")
    print(f"Relative overlap in CeCl6(2-): {orbital_overlap_CeCl6}")
    print("-" * 30)

    if orbital_overlap_CeF6 > orbital_overlap_CeCl6:
        conclusion = "stronger"
    else:
        # This case includes weaker or equal, but the problem states it's greater.
        conclusion = "weaker or equal"
    
    print(f"Conclusion: Given its higher orbital overlap, CeF6(2-) would display {conclusion} covalency compared to CeCl6(2-).")

if __name__ == '__main__':
    analyze_covalency()
