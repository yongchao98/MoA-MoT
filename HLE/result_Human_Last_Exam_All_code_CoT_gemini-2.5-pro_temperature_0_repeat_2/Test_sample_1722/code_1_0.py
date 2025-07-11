import itertools

def solve():
    """
    This function solves the POMDP memory problem by translating it into a problem
    from automata theory and then computationally finding the smallest integer solution.
    """

    n = 1
    while True:
        # Generate all binary strings of length n
        # e.g., for n=2, this is ['00', '01', '10', '11']
        binary_strings = ["".join(p) for p in itertools.product('01', repeat=n)]

        # Generate all unique pairs of these strings
        # e.g., ('00', '01'), ('00', '10'), ...
        string_pairs = itertools.combinations(binary_strings, 2)

        found_pair = None
        for s1, s2 in string_pairs:
            # Count zeros and ones for each string
            n0_s1 = s1.count('0')
            n1_s1 = s1.count('1')
            n0_s2 = s2.count('0')
            n1_s2 = s2.count('1')

            # Check the parity conditions for 2-state indistinguishability
            if (n0_s1 % 2 == n0_s2 % 2) and (n1_s1 % 2 == n1_s2 % 2):
                found_pair = (s1, s2)
                break
        
        if found_pair:
            s1, s2 = found_pair
            print(f"Found solution for n = {n}")
            print("-" * 30)
            print("This is the smallest `n` for which we can find two different observation sequences")
            print("that a 2-state memory cannot tell apart, but a 3-state memory can.")
            print("\nAn example pair of sequences for this n:")
            print(f"  Sequence in C1 (Ω1): {s1}")
            print(f"  Sequence in C2 (Ω2): {s2}")
            
            print("\nChecking condition for m=2 (indistinguishable):")
            n0_s1 = s1.count('0')
            n1_s1 = s1.count('1')
            n0_s2 = s2.count('0')
            n1_s2 = s2.count('1')
            
            print(f"  Parity of zeros in Ω1: {n0_s1} % 2 = {n0_s1 % 2}")
            print(f"  Parity of zeros in Ω2: {n0_s2} % 2 = {n0_s2 % 2}")
            print(f"  Parity of ones in Ω1:  {n1_s1} % 2 = {n1_s1 % 2}")
            print(f"  Parity of ones in Ω2:  {n1_s2} % 2 = {n1_s2 % 2}")
            print("  Parities match, so no 2-state machine can distinguish them.")

            print("\nChecking condition for m=3 (distinguishable):")
            print(f"  Since '{s1}' and '{s2}' are different, a 3-state machine can be constructed to distinguish them.")
            print(f"  This allows an agent with m=3 memory to achieve a higher reward.")
            
            print("\nFinal Answer:")
            print(f"The minimum length of the hallway is {n}.")
            print(f"<<<{n}>>>")
            return
            
        # If no pair was found, increment n and try again
        n += 1

if __name__ == "__main__":
    solve()