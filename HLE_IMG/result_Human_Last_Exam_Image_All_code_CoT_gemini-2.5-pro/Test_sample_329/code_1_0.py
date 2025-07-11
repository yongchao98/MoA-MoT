import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string stream
sys.stdout = captured_output = io.StringIO()


def solve_puzzle():
    """
    Analyzes four triple-slit simulations to find a four-digit code.

    The function deduces the behavior of three parameters (width 'w', wave number 'u', height 'h')
    for three slits in four different simulations based on visual cues from the provided images.
    It identifies the unique parameter for each simulation as per the problem's rules and
    constructs the final four-digit sequence from their corresponding codes.
    """
    # Define the encoding for each parameter and slit
    encoding = {
        'w1': 1, 'w2': 2, 'w3': 3,  # width (left, middle, right)
        'u1': 4, 'u2': 5, 'u3': 6,  # wave number (left, middle, right)
        'h1': 7, 'h2': 8, 'h3': 9,  # height (left, middle, right)
    }

    print("Step-by-step analysis of each simulation:")
    print("="*40)

    # --- Simulation 1 Analysis ---
    print("Analysis for Simulation 1:")
    print("- Brightness (h): The brightness at the slits increases from left to right. Therefore, 'h' is the monotonic parameter.")
    print("- Beam Spread (w): The middle beam is wider than the left and right beams, which look similar. A wider beam implies a narrower slit. Thus, the middle slit's width (w2) is unique. 'w' is the unique parameter.")
    print("- Steering (u): By elimination, 'u' must be the identical parameter, which is consistent with the symmetric, un-steered pattern.")
    unique_param_1 = 'w2'
    code_1 = encoding[unique_param_1]
    print(f"Conclusion: The unique parameter is w2 (middle slit width). Code: {code_1}")
    print("="*40)

    # --- Simulation 2 Analysis ---
    print("Analysis for Simulation 2:")
    print("- Brightness (h): The left slit is significantly brighter than the middle and right slits, which appear equally bright. Thus, the left slit's height (h1) is unique. 'h' is the unique parameter.")
    print("- Beam Spread (w): The beam width increases from left to right (narrowest on left, widest on right). This means slit width 'w' decreases monotonically from left to right.")
    print("- Steering (u): By elimination, 'u' is the identical parameter, consistent with the symmetric pattern.")
    unique_param_2 = 'h1'
    code_2 = encoding[unique_param_2]
    print(f"Conclusion: The unique parameter is h1 (left slit height). Code: {code_2}")
    print("="*40)

    # --- Simulation 3 Analysis ---
    print("Analysis for Simulation 3:")
    print("- Brightness (h): All three slits have the same initial brightness. 'h' is the identical parameter.")
    print("- Steering (u): The entire pattern is steered to the left, indicating a monotonic change in wave number 'u' across the slits.")
    print("- Beam Spread (w): The middle beam is narrower than the left and right beams, which look similar. A narrower beam implies a wider slit. Thus, the middle slit's width (w2) is unique.")
    unique_param_3 = 'w2'
    code_3 = encoding[unique_param_3]
    print(f"Conclusion: The unique parameter is w2 (middle slit width). Code: {code_3}")
    print("="*40)

    # --- Simulation 4 Analysis ---
    print("Analysis for Simulation 4:")
    print("- Beam Spread (w): All three beams have the same width. 'w' is the identical parameter.")
    print("- Brightness (h): The brightness at the slits decreases from left to right. 'h' is the monotonic parameter.")
    print("- Steering (u): By elimination, 'u' is the unique parameter. The pattern is steered to the right, and the left beam bends more than the others, indicating u1 is different.")
    unique_param_4 = 'u1'
    code_4 = encoding[unique_param_4]
    print(f"Conclusion: The unique parameter is u1 (left slit wave number). Code: {code_4}")
    print("="*40)

    # --- Final Answer ---
    final_sequence = f"{code_1}{code_2}{code_3}{code_4}"
    print("\nFinal Result Compilation:")
    print(f"The unique parameter codes are:")
    print(f"Simulation 1: {code_1}")
    print(f"Simulation 2: {code_2}")
    print(f"Simulation 3: {code_3}")
    print(f"Simulation 4: {code_4}")
    print(f"\nThe resulting four-digit sequence is {final_sequence}.")

solve_puzzle()
# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)
print("<<<2724>>>")