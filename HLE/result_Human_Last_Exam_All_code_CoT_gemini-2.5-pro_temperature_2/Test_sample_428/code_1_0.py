import sys

def solve():
    """
    This function calculates the expected time until the second collision of five particles.
    The user is expected to provide four positive integers N1, M1, N2, M2 as command-line arguments.
    """
    if len(sys.argv) != 5:
        print("Usage: python your_script_name.py N1 M1 N2 M2")
        print("Please provide four positive integer arguments for N1, M1, N2, M2.")
        return

    try:
        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])
        if not all(x > 0 for x in [N1, M1, N2, M2]):
            raise ValueError("All inputs must be positive integers.")
    except ValueError as e:
        print(f"Error: Invalid input. {e}")
        return

    # The initial gaps between adjacent particles are D1=N1, D2=M1, D3=N2, D4=M2.
    # The expected time E[tau] for the second collision is given by the formula:
    # E[tau] = (1/4) * sum_{1<=i<j<=4} D_i * D_j
    
    term12 = N1 * M1
    term13 = N1 * N2
    term14 = N1 * M2
    term23 = M1 * N2
    term24 = M1 * M2
    term34 = N2 * M2
    
    numerator = term12 + term13 + term14 + term23 + term24 + term34
    
    expected_time = numerator / 4.0
    
    print(f"The initial gaps are N1={N1}, M1={M1}, N2={N2}, M2={M2}.")
    print("The expectation of tau is calculated by the formula:")
    print(f"E[tau] = (N1*M1 + N1*N2 + N1*M2 + M1*N2 + M1*M2 + N2*M2) / 4")
    print(f"E[tau] = ({N1}*{M1} + {N1}*{N2} + {N1}*{M2} + {M1}*{N2} + {M1}*{M2} + {N2}*{M2}) / 4")
    print(f"E[tau] = ({term12} + {term13} + {term14} + {term23} + {term24} + {term34}) / 4")
    print(f"E[tau] = {numerator} / 4")
    print(f"E[tau] = {expected_time}")

if __name__ == "__main__":
    solve()
