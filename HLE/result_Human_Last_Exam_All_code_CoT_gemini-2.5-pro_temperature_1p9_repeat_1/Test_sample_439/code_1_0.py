def find_shortest_program_length(target_n):
    """
    This function simulates the computation of K(n) for a primitive recursive language.

    Since the language P is primitive recursive, every program is guaranteed to halt.
    This allows us to find K(n) with a simple brute-force search.

    This code simulates that search on a predefined, finite set of mock programs.
    """

    # In a real scenario, an interpreter for the language 'P' would generate
    # and run programs. Here, we mock a set of programs and their outputs.
    # The keys are program strings and values are their outputs.
    # Notice multiple programs can produce the same output (e.g., for 100).
    mock_p_programs = {
        'A': 0,
        'B': 1,
        'C': 5,
        'INIT': 10,
        'ADD1': 2,
        'LOOP': 100,
        'XYZ': 1000,
        'RECURSE': 2,
        'ABCDEFG': 5,
        'PROG_FINAL': 100,
        'A_REALLY_LONG_PROGRAM_NAME': 1000,
    }

    print(f"Searching for the shortest program that outputs the number {target_n}...")
    
    shortest_length = -1
    shortest_program = ""

    # We determine the maximum possible length from our mock set for the simulation.
    # In a real computation, this loop would be unbounded (for length in 1, 2, 3, ...).
    max_len = max(len(p) for p in mock_p_programs.keys()) if mock_p_programs else 0

    # 1. Iterate through all possible program lengths, from 1 upwards.
    for length in range(1, max_len + 1):
        print(f"\nChecking all programs of length {length}...")
        
        # 2. Find all programs of the current length.
        programs_of_current_length = {
            prog: output
            for prog, output in mock_p_programs.items()
            if len(prog) == length
        }

        if not programs_of_current_length:
            print("No programs of this length exist in our mock set.")
            continue

        # 3. "Execute" each program and check its output.
        for program, output in programs_of_current_length.items():
            print(f"  - Running program '{program}'... it outputs {output}.")
            if output == target_n:
                # 4. If a program produces the target, we have found the shortest length.
                shortest_length = length
                shortest_program = program
                print("\n-----------------------------------------------------")
                print(f"Success! Found a program that outputs {target_n}.")
                # Print the 'final equation' components as requested.
                print(f"Number (n): {target_n}")
                print(f"Shortest Length Found (K(n)): {shortest_length}")
                print(f"The final equation is: K({target_n}) = {shortest_length}")
                print("-----------------------------------------------------")
                return

    # This part would only be reached if no program could produce n.
    # However, for any n, a primitive recursive program can be constructed.
    print(f"\nCould not find any program that outputs {target_n} in our mock set.")

# Let's find K(100) using our simulation.
# Based on our mock set, 'LOOP' and 'PROG_FINAL' both output 100, but 'LOOP' has length 4,
# which is shorter than 'PROG_FINAL' (length 10). The search should find 'LOOP'.
find_shortest_program_length(100)