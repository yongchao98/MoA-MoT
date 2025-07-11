import sys
import io

# Helper to capture print output for clean final presentation
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_and_explain():
    """
    Analyzes the four triple-slit simulations to find the unique parameter
    for each and generates the corresponding four-digit code.
    """
    
    # The encoding map as specified in the problem description
    code_map = {
        'w1': 1, 'w2': 2, 'w3': 3,  # width of left, middle, right slit
        'u1': 4, 'u2': 5, 'u3': 6,  # wave number of left, middle, right slit
        'h1': 7, 'h2': 8, 'h3': 9   # height of left, middle, right slit
    }

    # This dictionary stores the result of the visual analysis for each simulation.
    # The logic is based on applying physics principles to the images.
    analysis_results = {
        1: {
            "unique_parameter": "h1",
            "reasoning": {
                "Unique (Rule 3)": "Parameter 'h' is unique. Slit 1 is dimmer than slits 2 and 3, which have similar brightness. Thus, h1 is the unique value.",
                "Monotonic (Rule 2)": "Parameter 'u' is monotonic. The entire pattern is steered to the left, which implies a monotonic phase change across the slits.",
                "Identical (Rule 1)": "Parameter 'w' is identical. The angular spread of the beams from all three slits appears to be the same."
            }
        },
        2: {
            "unique_parameter": "w2",
            "reasoning": {
                "Unique (Rule 3)": "Parameter 'w' is unique. The beam from the middle slit is wider than the other two, which implies slit 2 is narrower. Thus, w2 is the unique value.",
                "Identical (Rule 1)": "Parameter 'h' is identical. All three slits have the same initial brightness.",
                "Monotonic (Rule 2)": "Parameter 'u' is monotonic. The pattern is steered to the right."
            }
        },
        3: {
            "unique_parameter": "h2",
            "reasoning": {
                "Unique (Rule 3)": "Parameter 'h' is unique. The middle slit is dimmer than the two outer slits. Thus, h2 is the unique value.",
                "Identical (Rule 1)": "Parameter 'w' is identical. All three beams have similar angular spreads.",
                "Monotonic (Rule 2)": "Parameter 'u' is monotonic. The pattern is steered to the right."
            }
        },
        4: {
            "unique_parameter": "w3",
            "reasoning": {
                "Unique (Rule 3)": "Parameter 'w' is unique. The beam from the rightmost slit (3) is wider than the other two, implying slit 3 is narrower. Thus, w3 is the unique value.",
                "Identical (Rule 1)": "Parameter 'h' is identical. All three slits show the same initial brightness.",
                "Monotonic (Rule 2)": "Parameter 'u' is monotonic. The pattern is steered to the left."
            }
        }
    }

    final_code_digits = []

    print("Step-by-step analysis and code generation:\n")

    for i in range(1, 5):
        sim_data = analysis_results[i]
        unique_param = sim_data["unique_parameter"]
        code = code_map[unique_param]
        final_code_digits.append(str(code))

        print(f"--- Simulation {i} ---")
        print(f"Identified Unique Parameter: '{unique_param}'")
        print("Reasoning:")
        for rule, explanation in sim_data["reasoning"].items():
            print(f"  - {rule}: {explanation}")
        print(f"Code for '{unique_param}': {code}\n")

    # Assemble the final answer
    digit1, digit2, digit3, digit4 = final_code_digits
    
    print("="*30)
    print("Final Code Assembly")
    print("="*30)
    print(f"Simulation 1 Code: {digit1}")
    print(f"Simulation 2 Code: {digit2}")
    print(f"Simulation 3 Code: {digit3}")
    print(f"Simulation 4 Code: {digit4}")
    print(f"The final four-digit sequence is: {digit1}{digit2}{digit3}{digit4}")

solve_and_explain()
sys.stdout = old_stdout
print(captured_output.getvalue())