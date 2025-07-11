import numpy as np
import itertools

def check_diversity(code_name):
    """
    Checks the diversity order of a given space-time code by searching for
    singular difference matrices over a 4-QAM constellation.

    Args:
        code_name (str): The name of the code ('Sa', 'Sb', or 'Sc').

    Returns:
        int: The determined diversity order (1 or 2).
    """
    print(f"--- Analyzing Code {code_name} ---")
    
    # 4-QAM constellation (unit average energy)
    qam_symbols = [
        (1 + 1j) / np.sqrt(2),
        (1 - 1j) / np.sqrt(2),
        (-1 + 1j) / np.sqrt(2),
        (-1 - 1j) / np.sqrt(2)
    ]
    
    # Generate all possible pairs of input symbols (x1, x2)
    symbol_pairs = list(itertools.product(qam_symbols, repeat=2))
    
    # Iterate through all pairs of distinct codewords
    # A codeword is defined by a symbol_pair (x1, x2)
    for s1_symbols in symbol_pairs:
        for s2_symbols in symbol_pairs:
            if s1_symbols == s2_symbols:
                continue

            x1, x2 = s1_symbols
            x1_prime, x2_prime = s2_symbols

            delta_x1 = x1 - x1_prime
            delta_x2 = x2 - x2_prime

            delta_S = np.zeros((2, 2), dtype=np.complex128)
            det = 0

            if code_name == 'Sa':
                # S_a = [[x1, x2], [x2, x1]]
                delta_S = np.array([[delta_x1, delta_x2],
                                    [delta_x2, delta_x1]])
                det = delta_x1**2 - delta_x2**2
            elif code_name == 'Sb':
                # S_b = [[x1, x2], [x2, x1*]]
                delta_S = np.array([[delta_x1, delta_x2],
                                    [delta_x2, np.conj(delta_x1)]])
                det = delta_x1 * np.conj(delta_x1) - delta_x2**2
            elif code_name == 'Sc':
                # S_c = [[-x1*, x2], [-x2*, -x1]]
                delta_S = np.array([[-np.conj(delta_x1), delta_x2],
                                    [-np.conj(delta_x2), -delta_x1]])
                det = (-np.conj(delta_x1)) * (-delta_x1) - delta_x2 * (-np.conj(delta_x2))
                det = np.abs(delta_x1)**2 + np.abs(delta_x2)**2
            
            # Check if the determinant is zero
            if np.isclose(det, 0):
                print(f"Found a singular difference matrix for code {code_name}!")
                print(f"  x1 = {x1:.3f}, x2 = {x2:.3f}")
                print(f"  x1' = {x1_prime:.3f}, x2' = {x2_prime:.3f}")
                print(f"  Delta_x1 = {delta_x1:.3f}, Delta_x2 = {delta_x2:.3f}")
                if code_name == 'Sa':
                    print(f"  det(Delta_S) = (Delta_x1)^2 - (Delta_x2)^2 = ({delta_x1:.3f})^2 - ({delta_x2:.3f})^2 = {det:.3f}")
                elif code_name == 'Sb':
                    print(f"  det(Delta_S) = |Delta_x1|^2 - (Delta_x2)^2 = |{delta_x1:.3f}|^2 - ({delta_x2:.3f})^2 = {np.abs(delta_x1)**2:.3f} - ({delta_x2**2:.3f}) = {det:.3f}")
                print("Conclusion: Diversity order is 1.\n")
                return 1

    print(f"No singular difference matrices found for code {code_name}.")
    print("Conclusion: The code has full diversity. Diversity order is 2.\n")
    return 2

if __name__ == '__main__':
    diversity_orders = {}
    
    # (a) Calculate diversity order for each code
    diversity_orders['Sa'] = check_diversity('Sa')
    diversity_orders['Sb'] = check_diversity('Sb')
    diversity_orders['Sc'] = check_diversity('Sc')
    
    print("--- Summary of Results ---")
    print("(a) Diversity Orders:")
    for code, order in diversity_orders.items():
        print(f"  - Diversity order for code {code} is: {order}")
    
    # (b) Find the code with the maximum directivity (highest diversity)
    max_diversity = 0
    best_code = None
    for code, order in diversity_orders.items():
        if order > max_diversity:
            max_diversity = order
            best_code = code

    print("\n(b) Maximum Directivity:")
    print(f"Code {best_code} has the highest diversity order ({max_diversity}).")
    print("Therefore, it provides the best performance, interpreted as maximum directivity.")
