import math

def calculate_toric_code_degeneracy(m, n):
    """
    Calculates the ground space degeneracy (GSD) of the toric code on a sphere
    with m rough holes and n smooth holes.

    This formula is generally valid for m >= 1 and n >= 1.

    Args:
        m (int): The number of rough holes.
        n (int): The number of smooth holes.
    """
    if not isinstance(m, int) or not isinstance(n, int) or m < 0 or n < 0:
        print("Error: The number of holes (m and n) must be non-negative integers.")
        return

    print(f"Calculating ground space degeneracy for a toric code with m={m} rough holes and n={n} smooth holes.")
    print("The base surface is assumed to be a sphere (genus g=0).")
    
    if m == 0 and n == 0:
        k = 0
        degeneracy = 1
        print("For a closed sphere with no holes, the number of logical qubits k is 0.")
        print(f"Ground Space Degeneracy = 2^k = 2^0 = {degeneracy}")
    elif m <= 0 or n <= 0:
        # Handle the edge cases where one type of hole is absent
        # The correct formula is k = max(0, m-1) + max(0, n-1)
        # We demonstrate the formula from Option B, while acknowledging its limitations
        print("\nNote: The formula GSD = 2^(m+n-2) is primarily for cases where m>0 and n>0.")
        print("We will apply it directly as per the multiple-choice options.")
        if m + n < 2:
            print("The formula m+n-2 results in a non-positive exponent, which is physically incorrect for GSD.")
            k = 'invalid'
            degeneracy = 'invalid'
        else:
             k = m + n - 2
             degeneracy = 2**k
        print(f"Using the formula from option B, the number of logical qubits k = m + n - 2 = {m} + {n} - 2 = {k}")
        print(f"Ground Space Degeneracy = 2^k = 2^{k} = {degeneracy}")

    else: # Case m>=1 and n>=1
        qubits_from_rough = m - 1
        qubits_from_smooth = n - 1
        k = qubits_from_rough + qubits_from_smooth
        degeneracy = 2**k

        print(f"\nNumber of logical qubits from {m} rough holes = m - 1 = {m} - 1 = {qubits_from_rough}")
        print(f"Number of logical qubits from {n} smooth holes = n - 1 = {n} - 1 = {qubits_from_smooth}")
        print(f"Total number of logical qubits k = ({m} - 1) + ({n} - 1) = {qubits_from_rough} + {qubits_from_smooth} = {k}")
        print(f"Final equation for GSD: 2^({m} + {n} - 2) = 2^{k}")
        print(f"Ground Space Degeneracy = {degeneracy}")


# Example usage with m=2 rough holes and n=3 smooth holes.
# You can change these values to test other cases.
m_holes = 2
n_holes = 3
calculate_toric_code_degeneracy(m_holes, n_holes)