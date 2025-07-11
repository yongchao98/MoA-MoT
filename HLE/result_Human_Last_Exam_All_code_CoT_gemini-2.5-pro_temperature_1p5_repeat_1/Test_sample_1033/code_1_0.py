import sys

def solve():
    """
    This function decodes the pattern in the sequence and calculates the next three terms.
    The relationship between the letters L1, L2, L3 in each triplet is defined by the formula:
    L3 = (a*L1 + b*L2 + c*i + d) % 26
    where i is the 0-based index of the triplet within its group.
    """

    def to_int(c):
        """Converts a character 'A'-'Z' to an integer 0-25."""
        return ord(c) - ord('A')

    def to_char(n):
        """Converts an integer 0-25 to a character 'A'-'Z'."""
        return chr(n + ord('A'))

    # Based on analysis of the sequence, the coefficients are determined to be:
    a, b, c, d = 13, 7, 14, 13
    
    # The last group in the sequence starts with 'N', so the next group starts with 'O'.
    l1_char = 'O'
    l1_val = to_int(l1_char)
    
    # We predict the L2 sequence for the 'O' group starts with 'A', 'B', 'C'.
    next_l2_chars = ['A', 'B', 'C']
    
    results = []

    print("Calculating the next three terms based on the derived formula:")
    print(f"Formula: L3 = ({a}*L1 + {b}*L2 + {c}*i + {d}) % 26\n")

    # Calculate the next three terms
    for i, l2_char in enumerate(next_l2_chars):
        l2_val = to_int(l2_char)
        
        # Apply the formula
        l3_val = (a * l1_val + b * l2_val + c * i + d) % 26
        l3_char = to_char(l3_val)
        
        term = f"{l1_char}{l2_char}{l3_char}"
        results.append(term)
        
        # Print the calculation for transparency
        print(f"For i={i}, L1='{l1_char}'({l1_val}), L2='{l2_char}'({l2_val}):")
        print(f"  L3 = ({a}*{l1_val} + {b}*{l2_val} + {c}*{i} + {d}) % 26 = {l3_val} ('{l3_char}')")
        print(f"  --> Resulting term: {term}\n")
    
    final_answer = " ".join(results)
    return final_answer

# Execute the solution
final_answer = solve()

# The final result in the requested format
sys.stdout.write(f'<<<{final_answer}>>>\n')