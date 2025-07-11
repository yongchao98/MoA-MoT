def solve():
    """
    This function calculates the next 10 characters of the sequence based on the
    identified rule: Thue-Morse sequence applied to the Moser-de Bruijn sequence.
    """
    prefix_len = 4671
    results = []
    
    start_n = prefix_len + 1
    end_n = prefix_len + 10

    for n in range(start_n, end_n + 1):
        # Step 1: Calculate the Moser-de-Bruijn number for n (1-indexed).
        # This is equivalent to interpreting the binary representation of n
        # as a number in base 4.
        n_binary_str = bin(n)[2:]
        moser_de_bruijn_num = int(n_binary_str, 4)

        # Step 2: Calculate the Thue-Morse value of the resulting number.
        # This is the parity of the count of '1's in its binary representation.
        pop_count = bin(moser_de_bruijn_num).count('1')
        thue_morse_val = pop_count % 2
        
        results.append(str(thue_morse_val))
        
    print("".join(results))

solve()