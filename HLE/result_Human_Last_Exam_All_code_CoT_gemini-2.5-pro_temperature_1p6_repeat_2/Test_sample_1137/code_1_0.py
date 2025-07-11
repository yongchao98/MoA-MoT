import sys
import io

# Buffer to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_litigation_forum():
    """
    Analyzes the litigation options for the RE1 vs. RE2 dispute and identifies the best forum.
    """
    # Key requirements derived from the problem description
    # 1. Is it an originating court (not an appellate court)?
    # 2. Can it handle large monetary values (>$35k)?
    # 3. Is it specialized for complex commercial matters?
    # 4. Is it designed for a speedy resolution?
    # 5. Is it the correct jurisdiction (provincial for this case)?

    options = {
        'A': {'name': 'Ontario Court of Appeal', 'originating': 0, 'large_value': 1, 'commercial_spec': 1, 'speedy': 0, 'correct_jurisdiction': 1},
        'B': {'name': 'Commercial List', 'originating': 1, 'large_value': 1, 'commercial_spec': 1, 'speedy': 1, 'correct_jurisdiction': 1},
        'C': {'name': 'Superior Court of Justice', 'originating': 1, 'large_value': 1, 'commercial_spec': 0, 'speedy': 0, 'correct_jurisdiction': 1},
        'D': {'name': 'Small Claims Court', 'originating': 1, 'large_value': 0, 'commercial_spec': 0, 'speedy': 1, 'correct_jurisdiction': 1},
        'E': {'name': 'Federal Court of Canada', 'originating': 1, 'large_value': 1, 'commercial_spec': 0, 'speedy': 0, 'correct_jurisdiction': 0}
    }

    best_option_key = None
    max_score = -1
    best_option_details = {}

    print("Evaluating litigation forum options based on suitability (1=suitable, 0=unsuitable):\n")

    for key, details in options.items():
        score = details['originating'] + details['large_value'] + details['commercial_spec'] + details['speedy'] + details['correct_jurisdiction']
        print(f"Option {key} ({details['name']}): Score = {score}")
        if score > max_score:
            max_score = score
            best_option_key = key
            best_option_details = details

    print("\n--- Analysis Conclusion ---")
    print(f"The best option is '{best_option_key}' with a score of {max_score}.")
    print("This is because the dispute is a complex, high-value commercial matter where speed is a key priority, making the specialized Commercial List the most appropriate forum.")
    
    # Fulfilling the request to show the final equation with numbers
    print("\nFinal equation for the best option's score:")
    equation_parts = [
        f"{best_option_details['originating']} (Originating Court)",
        f"{best_option_details['large_value']} (Large Value)",
        f"{best_option_details['commercial_spec']} (Commercial Specialty)",
        f"{best_option_details['speedy']} (Speedy Process)",
        f"{best_option_details['correct_jurisdiction']} (Correct Jurisdiction)"
    ]
    print(" + ".join(equation_parts) + f" = {max_score}")

    print("\n<<<B>>>")

# Run the analysis
solve_litigation_forum()

# Get the captured output and print it to the actual console
sys.stdout = old_stdout
print(captured_output.getvalue())