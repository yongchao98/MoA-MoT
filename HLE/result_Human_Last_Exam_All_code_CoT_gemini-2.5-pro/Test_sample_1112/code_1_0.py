import sys

def generate_sequence(bit_limit):
    """
    Generates the Chaitin-Tromp sequence (OEIS A102235).

    This algorithm simulates a universal prefix-free Turing machine.
    'progs' is the set of non-halted programs.
    'p' and 'q' are a pair of prefixes from a complete prefix code.
    The loop decides whether to extend a program (if 'p' wins) or
    to halt one and output a bit (if 'q' wins).
    """
    # The set of "running" programs, which form a prefix-free code.
    progs = {''}
    # 'n' is used to generate the next pair of prefixes 'p' and 'q'.
    n = 1
    output_bits = []

    # 'z' is used as a placeholder that is lexicographically larger than any binary string.
    # This simplifies finding the minimum program prefixed by p or q.
    Z_SENTINEL = 'z'

    while len(output_bits) < bit_limit:
        # Generate the next prefix pair (p, q) from integer n.
        # e.g., n=1 -> p='1', q='0'; n=2 -> p='10', q='11'
        p = bin(n)[2:]
        q = p[:-1] + ('1' if p[-1] == '0' else '0')

        # Find the lexicographically smallest program in our set that starts with p.
        p_prefixed = [s for s in progs if s.startswith(p)]
        p_choice = min(p_prefixed) if p_prefixed else Z_SENTINEL

        # Find the lexicographically smallest program in our set that starts with q.
        q_prefixed = [s for s in progs if s.startswith(q)]
        q_choice = min(q_prefixed) if q_prefixed else Z_SENTINEL

        # The 'winner' is the one with the lexicographically smaller choice.
        if p_choice < q_choice:
            # p wins: Extend the program 'p_choice' by adding '0' and '1'.
            # This maintains the prefix-free property of the set.
            progs.remove(p_choice)
            progs.add(p_choice + '0')
            progs.add(p_choice + '1')
            n += 1  # Move to the next pair of prefixes.
        else:
            # q wins: The program 'q_choice' halts. We output a bit.
            # The output bit is the character in q_choice immediately after the prefix q.
            bit = q_choice[len(q)]
            output_bits.append(bit)
            
            # The program has halted, so we remove it from the set.
            progs.remove(q_choice)
            # We do *not* increment n, as the same (p,q) pair must be tested
            # against the newly modified set of programs.

    return "".join(output_bits)

# The user's prefix is 4096 characters long. We need the next 10.
prefix_length = 4096
num_next_chars = 10
total_length = prefix_length + num_next_chars

# Generate the sequence.
# Note: This is computationally intensive and may take a few seconds.
full_sequence = generate_sequence(total_length)

# Extract the 10 characters that follow the prefix of length 4096.
next_10_chars = full_sequence[prefix_length:]

print(next_10_chars)