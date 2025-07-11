import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new output stream
sys.stdout = new_stdout = io.StringIO()


def solve_mosquito_problem():
    """
    This function analyzes pond characteristics to determine which poses the greatest
    medical threat based on mosquito abundance.
    It assigns scores to size and age to quantify the threat level.
    """

    ponds = [
        {"id": "A", "description": "10 feet square, one year old", "size": 10, "age": 1},
        {"id": "C", "description": "30 feet square, one year old", "size": 30, "age": 1},
        {"id": "D", "description": "10 feet square, five years old", "size": 10, "age": 5},
        {"id": "E", "description": "30 feet square, five years old", "size": 30, "age": 5},
    ]

    best_pond = None
    max_threat_score = -1

    final_equation_components = {}

    print("Analyzing ponds based on ecological factors...\n")

    for pond in ponds:
        # Assign a score based on size (larger is better)
        # Larger area supports a larger population.
        size_score = 1 if pond["size"] == 10 else 2

        # Assign a score based on age (older is better)
        # Older ponds have more established ecosystems and food sources.
        age_score = 1 if pond["age"] == 1 else 2

        # Threat score is multiplicative, as the factors amplify each other.
        threat_score = size_score * age_score

        if threat_score > max_threat_score:
            max_threat_score = threat_score
            best_pond = pond
            final_equation_components = {
                'size_score': size_score,
                'age_score': age_score,
                'threat_score': threat_score
            }

    print("Conclusion: The greatest threat is posed by the pond with the largest size and oldest age.")
    print(f"The best option is Pond {best_pond['id']}: {best_pond['description']}.")
    print("\nCalculation for the highest threat score:")
    # Final step to show the numbers in the final equation
    size_num = final_equation_components['size_score']
    age_num = final_equation_components['age_score']
    threat_num = final_equation_components['threat_score']

    print(f"Size Score ({best_pond['size']} ft) = {size_num}")
    print(f"Age Score ({best_pond['age']} years) = {age_num}")
    print(f"Final Threat Score = {size_num} (Size Score) * {age_num} (Age Score) = {threat_num}")

solve_mosquito_problem()

# Get the content from the new_stdout
output = new_stdout.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the captured output
print(output)
print("<<<E>>>")