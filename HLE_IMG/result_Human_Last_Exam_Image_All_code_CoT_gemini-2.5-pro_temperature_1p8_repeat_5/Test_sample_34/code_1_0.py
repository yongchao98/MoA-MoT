import sys
import io

# A helper function to print explanations in a structured way.
def explain(text):
    print(text)

def main():
    """
    Interprets a Lévi-Strauss kinship diagram and identifies the systems it represents.
    """
    # Step 1: Interpret the provided Lévi-Strauss kinship diagram.
    explain("Step 1: Interpreting the Lévi-Strauss Kinship Diagram")
    explain("="*50)
    explain("The diagram represents an 'atom of kinship' with the following relationships and attitudes:")
    explain("  - Symbols: Δ (male), o (female), = (marriage)")
    explain("  - Attitudes: + (positive/familiar), - (negative/formal)")
    explain("")

    # Based on the diagram:
    # 1. Brother(Δ) - Sister(o) -> The line between them is marked '-'
    # 2. Husband(right Δ) - Wife(o) -> The line between them is marked '+'
    # 3. Father(right Δ) - Son(bottom Δ) -> The line between them is marked '-'
    # 4. Mother's Brother(left Δ) - Sister's Son(bottom Δ) -> The line between them is marked '+'
    diagram_pattern = {
        "Husband-Wife": "+",
        "Brother-Sister": "-",
        "Father-Son": "-",
        "Mother's_Brother-Sister's_Son": "+"
    }
    
    explain("From the diagram, we deduce the following pattern of relationships:")
    print(f"  - Husband-Wife: {diagram_pattern['Husband-Wife']} (Familiar)")
    print(f"  - Brother-Sister: {diagram_pattern['Brother-Sister']} (Formal)")
    print(f"  - Father-Son: {diagram_pattern['Father-Son']} (Formal)")
    print(f"  - Mother’s Brother-Sister’s Son: {diagram_pattern['Mother\'s_Brother-Sister\'s_Son']} (Familiar)")
    explain("")

    # Step 2: Define the kinship patterns for the societies listed in the answer choices.
    explain("Step 2: Defining Known Kinship Patterns for Comparison")
    explain("="*50)
    society_patterns = {
        "Trobriand": {"lineage": "matrilineal", "pattern": {
            "Husband-Wife": "+", "Brother-Sister": "-", "Father-Son": "+", "Mother's_Brother-Sister's_Son": "-"}},
        "Siuoi": {"lineage": "matrilineal", "pattern": {
            "Husband-Wife": "+", "Brother-Sister": "-", "Father-Son": "+", "Mother's_Brother-Sister's_Son": "-"}},
        "Lake Kubutu": {"lineage": "patrilineal", "pattern": {
            "Husband-Wife": "+", "Brother-Sister": "-", "Father-Son": "-", "Mother's_Brother-Sister's_Son": "+"}},
        "Tonga": {"lineage": "patrilineal", "pattern": {
            "Husband-Wife": "+", "Brother-Sister": "-", "Father-Son": "-", "Mother's_Brother-Sister's_Son": "+"}},
        "Cherkess": {"lineage": "patrilineal", "pattern": {
            "Husband-Wife": "-", "Brother-Sister": "+", "Father-Son": "-", "Mother's_Brother-Sister's_Son": "+"}}
    }
    
    # Per instructions to output "each number in the final equation," we will display the values for each comparison.
    explain("Known patterns from anthropological studies:")
    for society, data in society_patterns.items():
        p = data['pattern']
        print(f"- {society} ({data['lineage']}): "
              f"HW={p['Husband-Wife']}, BS={p['Brother-Sister']}, FS={p['Father-Son']}, "
              f"MBSS={p['Mother\'s_Brother-Sister\'s_Son']}")
    explain("")

    # Step 3: Compare the diagram's pattern with the societies in each answer choice.
    explain("Step 3: Evaluating Answer Choices")
    explain("="*50)
    answer_choices = {
        "A": ["Trobriand", "Siuoi"],
        "B": ["Siuoi", "Lake Kubutu"],
        "C": ["Lake Kubutu", "Tonga"],
        "D": ["Tonga", "Cherkess"],
        "E": ["Cherkess", "Trobriand"]
    }

    correct_answer = ""
    for choice, societies in answer_choices.items():
        explain(f"Checking Choice {choice}: {societies[0]} and {societies[1]}")
        
        society1, society2 = societies[0], societies[1]
        match1 = society_patterns[society1]["pattern"] == diagram_pattern
        match2 = society_patterns[society2]["pattern"] == diagram_pattern
        
        print(f"  - Comparing diagram to {society1}: Match = {match1}")
        print(f"  - Comparing diagram to {society2}: Match = {match2}")

        if match1 and match2:
            correct_answer = choice
            explain(f"  Result: Both societies match. Choice {choice} is correct.\n")
        else:
            explain(f"  Result: This choice is incorrect.\n")

    explain("Final Conclusion:")
    explain(f"The analysis shows that the systems correctly represented by the diagram are {answer_choices[correct_answer][0]} and {answer_choices[correct_answer][1]}.")

if __name__ == "__main__":
    main()