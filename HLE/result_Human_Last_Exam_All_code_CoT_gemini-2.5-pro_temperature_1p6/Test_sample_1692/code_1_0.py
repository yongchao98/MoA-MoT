import sys

def find_intermediate_step():
    """
    Finds the intermediate step in a 1D elementary cellular automaton
    by testing all 256 possible rules.
    """
    gen_t0_str = "01101001"
    gen_t2_str = "10000111"

    # Convert string states to lists of integers
    gen_t0 = [int(c) for c in gen_t0_str]
    gen_t2 = [int(c) for c in gen_t2_str]
    n = len(gen_t0)
    
    found_solution = None

    for rule_num in range(256):
        # Convert the rule number into its 8-bit binary representation.
        # This string acts as a lookup table for the 8 possible neighborhoods.
        # e.g., rule 30 is '00011110'. The neighborhood '111' maps to the first
        # bit, '110' to the second, and so on.
        rule_bits = format(rule_num, '08b')

        # === Step 1: Calculate the intermediate state (gen_t1) from gen_t0 ===
        gen_t1_candidate = [0] * n
        for i in range(n):
            # Get neighborhood with periodic (wrapping) boundaries
            left_val = gen_t0[(i - 1 + n) % n]
            center_val = gen_t0[i]
            right_val = gen_t0[(i + 1) % n]
            
            # The neighborhood (e.g., (1,1,1)) forms a 3-bit number.
            # Its decimal value (e.g., 7) corresponds to an index in the rule_bits.
            # The standard mapping is (111->0, 110->1, ..., 000->7).
            # So, index = 7 - (left*4 + center*2 + right*1).
            neighborhood_val = 4 * left_val + 2 * center_val + 1 * right_val
            rule_index = 7 - neighborhood_val
            
            gen_t1_candidate[i] = int(rule_bits[rule_index])

        # === Step 2: Calculate the final state from our candidate gen_t1 ===
        gen_t2_calculated = [0] * n
        for i in range(n):
            # Get neighborhood from the candidate state
            left_val = gen_t1_candidate[(i - 1 + n) % n]
            center_val = gen_t1_candidate[i]
            right_val = gen_t1_candidate[(i + 1) % n]

            neighborhood_val = 4 * left_val + 2 * center_val + 1 * right_val
            rule_index = 7 - neighborhood_val
            
            gen_t2_calculated[i] = int(rule_bits[rule_index])

        # === Step 3: Check if the calculated state matches the known final state ===
        if gen_t2_calculated == gen_t2:
            # As the problem states there's only one valid solution,
            # we can stop once we find it.
            found_solution = "".join(map(str, gen_t1_candidate))
            break
            
    if found_solution:
        print(found_solution)
        # The line below is for the platform's specific answer format.
        sys.stdout.write(f'<<<__{found_solution}__>>>')

# Run the solver
find_intermediate_step()