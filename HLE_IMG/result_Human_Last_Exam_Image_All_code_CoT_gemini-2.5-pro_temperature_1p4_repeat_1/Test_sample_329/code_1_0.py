def solve_simulation_puzzle():
    """
    This function analyzes the four triple-slit simulations to find the unique parameter
    for each, according to the provided rules, and prints the final four-digit code.
    """

    # Encoding map for parameters: w(1-3), u(4-6), h(7-9)
    # w1->1, w2->2, w3->3, u1->4, u2->5, u3->6, h1->7, h2->8, h3->9
    encoding = {
        "w1": 1, "w2": 2, "w3": 3,
        "u1": 4, "u2": 5, "u3": 6,
        "h1": 7, "h2": 8, "h3": 9
    }

    # Analysis results for each simulation based on visual inspection
    # Each tuple contains: (parameter for Rule 1, parameter for Rule 2, parameter for Rule 3)
    analysis_results = {
        1: {'identical': 'w', 'monotonic': 'h', 'unique_param': 'u2'},
        2: {'identical': 'h', 'monotonic': 'w', 'unique_param': 'u1'},
        3: {'identical': 'w', 'monotonic': 'u', 'unique_param': 'h2'},
        4: {'identical': 'w', 'monotonic': 'u', 'unique_param': 'h2'}
    }

    final_sequence = []

    print("Step-by-step analysis:")
    for i in range(1, 5):
        unique_param_key = analysis_results[i]['unique_param']
        param_code = encoding[unique_param_key]
        final_sequence.append(str(param_code))
        
        print(f"\n--- Simulation {i} ---")
        if i == 1:
            print("Observation: Brightness increases left-to-right (h is monotonic). Beam widths are equal (w is identical).")
            print("Conclusion: The wave number (u) must follow Rule 3. The symmetric pattern implies the middle slit is unique.")
            print(f"Unique Parameter: {unique_param_key} -> Code: {param_code}")
        elif i == 2:
            print("Observation: Brightness is equal (h is identical). Beam width increases left-to-right, so slit width w is monotonic.")
            print("Conclusion: The wave number (u) must follow Rule 3. The rightward steer implies the left slit is unique.")
            print(f"Unique Parameter: {unique_param_key} -> Code: {param_code}")
        elif i == 3:
            print("Observation: The pattern is steered left (u is monotonic). Beam widths are equal (w is identical).")
            print("Conclusion: Height (h) must follow Rule 3. The middle slit is dimmer than the equal outer slits.")
            print(f"Unique Parameter: {unique_param_key} -> Code: {param_code}")
        elif i == 4:
            print("Observation: The pattern is steered left (u is monotonic). Beam widths are equal (w is identical).")
            print("Conclusion: Height (h) must follow Rule 3. The middle slit is brighter than the equal outer slits.")
            print(f"Unique Parameter: {unique_param_key} -> Code: {param_code}")

    # The problem asks for the four-digit sequence.
    final_code = "".join(final_sequence)
    print("\n-----------------------------------------")
    print(f"The unique parameter codes for simulations 1, 2, 3, and 4 are {final_sequence[0]}, {final_sequence[1]}, {final_sequence[2]}, and {final_sequence[3]}.")
    print(f"The final four-digit sequence is: {final_code}")


solve_simulation_puzzle()
<<<5488>>>