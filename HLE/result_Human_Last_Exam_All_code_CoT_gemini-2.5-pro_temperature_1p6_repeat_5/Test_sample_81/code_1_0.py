import sys
import io

# Redirect stdout to capture print output for the final result
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_contract_terms():
    """
    This function analyzes the sub-clauses of answer choice E.
    It assigns a "boilerplate score" (1-10) to each term, where 10 is
    a very common and expected term and 1 is a very surprising and
    unexpected term.
    The goal is to programmatically identify a "hidden" material term.
    """
    print("Analyzing terms from choice E for 'Reasonable Expectation' violation...\n")

    # A dictionary mapping clauses from choice E to a "boilerplate score"
    clauses = {
        "(i) no commercial use": 9,
        "(ii) no selling data": 9,
        "(iii) no copyright infringement": 10,
        "(iv) authorized users only": 9,
        "(v) no scraping": 10,
        "(vi) no use for Illinois residents": 1,
        "(vii) no bypassing security": 8
    }

    total_score = 0
    num_clauses = len(clauses)
    hidden_term_text = ""
    hidden_term_score = 11  # Initialize higher than max score

    print("Clause 'Boilerplate' Score Analysis:")
    for clause, score in clauses.items():
        print(f"Score: {score}/10 | Term: {clause}")
        total_score += score
        # Identify the most surprising term (the one with the lowest score)
        if score < hidden_term_score:
            hidden_term_score = score
            hidden_term_text = clause

    average_score = total_score / num_clauses
    
    print("\n--- Calculation ---")
    print("Equation for Average 'Boilerplate' Score:")
    # Demonstrating the final equation with each number
    score_list = [str(s) for s in clauses.values()]
    equation_str = f"({ ' + '.join(score_list) }) / {num_clauses}"
    print(f"{equation_str} = {total_score} / {num_clauses} = {average_score:.2f}\n")


    print("--- Conclusion ---")
    print(f"The most surprising term is '{hidden_term_text}' with a score of {hidden_term_score}.")
    print("Burying a significant service limitation (a ban on use for an entire state)")
    print("within a long list of standard terms is a classic example of hiding a material term.")
    print("This violates the 'doctrine of reasonable expectations', making E the most likely answer.")

analyze_contract_terms()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

print(output)
print("<<<E>>>")