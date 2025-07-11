import numpy as np

def analyze_stc():
    """
    Analyzes three space-time codes to determine their diversity order
    by examining the determinant of the codeword difference matrix.
    """
    # Use 4-QAM constellation symbols for demonstration
    qam_symbols = [1+1j, 1-1j, -1+1j, -1-1j]
    
    # Generate all possible codeword symbol pairs (x1, x2)
    symbol_pairs = []
    for x1 in qam_symbols:
        for x2 in qam_symbols:
            symbol_pairs.append((x1, x2))

    # Define the code matrix functions
    def get_Sa(symbols):
        x1, x2 = symbols
        return np.array([[x1, x2], [x2, x1]])

    def get_Sb(symbols):
        x1, x2 = symbols
        return np.array([[x1, x2], [x2, np.conj(x1)]])

    def get_Sc(symbols):
        x1, x2 = symbols
        return np.array([[-np.conj(x1), x2], [-np.conj(x2), -x1]])

    codes = {
        "S_a": (get_Sa, "det = e1^2 - e2^2"),
        "S_b": (get_Sb, "det = |e1|^2 - e2^2"),
        "S_c": (get_Sc, "det = |e1|^2 + |e2|^2")
    }

    diversity_orders = {}

    for name, (code_func, det_formula) in codes.items():
        print(f"--- Analyzing Code {name} ({det_formula}) ---")
        min_rank = 2  # Assume full rank initially
        found_rank_deficient = False
        
        num_pairs = len(symbol_pairs)
        for i in range(num_pairs):
            for j in range(num_pairs):
                if i == j:
                    continue  # Codewords must be distinct

                symbols1 = symbol_pairs[i]
                symbols2 = symbol_pairs[j]
                
                S1 = code_func(symbols1)
                S2 = code_func(symbols2)
                
                delta_S = S1 - S2
                det_delta_S = np.linalg.det(delta_S)

                # If determinant is close to zero, matrix is rank-deficient
                if np.isclose(det_delta_S, 0):
                    min_rank = 1
                    e1 = symbols1[0] - symbols2[0]
                    e2 = symbols1[1] - symbols2[1]
                    print("Found a rank-deficient case (Rank = 1).")
                    print(f"  Symbols 1 (x1, x2): ({symbols1[0]}, {symbols1[1]})")
                    print(f"  Symbols 2 (x1', x2'): ({symbols2[0]}, {symbols2[1]})")
                    print(f"  Error symbols (e1, e2): ({e1}, {e2})")
                    if name == "S_a":
                        print(f"  Checking determinant equation: e1^2 - e2^2 = ({e1})^2 - ({e2})^2 = {e1**2 - e2**2:.1f}")
                    elif name == "S_b":
                        print(f"  Checking determinant equation: |e1|^2 - e2^2 = |{e1}|^2 - ({e2})^2 = {abs(e1)**2 - e2**2:.1f}")
                    found_rank_deficient = True
                    break
            if found_rank_deficient:
                break
        
        diversity_orders[name] = min_rank
        if not found_rank_deficient:
            print("No rank-deficient cases found. The code has full rank (Rank = 2).")
            # Show an example calculation for S_c
            symbols1 = (1+1j, 1-1j)
            symbols2 = (-1+1j, -1-1j)
            e1 = symbols1[0] - symbols2[0]
            e2 = symbols1[1] - symbols2[1]
            print("Example calculation for a distinct pair:")
            print(f"  e1={e1}, e2={e2}")
            print(f"  det = |e1|^2 + |e2|^2 = |{e1}|^2 + |{e2}|^2 = {abs(e1)**2:.1f} + {abs(e2)**2:.1f} = {abs(e1)**2 + abs(e2)**2:.1f}")
        
        print(f"Result: The diversity order for {name} is {min_rank}.\n")

    print("--- Final Conclusion ---")
    print("(a) What is the diversity order for each code?")
    for name, order in diversity_orders.items():
        print(f"    - Code {name}: {order}")
    
    print("\n(b) Which code provides the maximum directivity?")
    print("Code S_c provides the maximum directivity by achieving the full diversity order of 2.")

analyze_stc()