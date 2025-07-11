import math

def calculate_convolution_operations():
    """
    Calculates the number of DFT+IDFT operations for overlap-add and
    overlap-save methods for linear convolution.
    """
    # Step 1: Define Parameters
    L_x = 1200  # Length of the long sequence x[n]
    M = 90      # Length of the short sequence h[n]
    N = 128     # DFT size

    print("Step 1: Define Given Parameters")
    print(f"  - Length of sequence x[n] (L_x): {L_x}")
    print(f"  - Length of sequence h[n] (M): {M}")
    print(f"  - DFT size (N): {N}")
    print("-" * 50)

    # Step 2: Calculate the Effective Data Block Size (L)
    # The condition to avoid aliasing is N >= L + M - 1.
    # We choose the largest L, so L = N - M + 1.
    # This L is the number of new/valid samples per block.
    L = N - M + 1

    print("Step 2: Calculate Effective Data Block Size (L)")
    if L <= 0:
        print("Error: DFT size N is too small for the given filter length M.")
        print("N must be greater than M-1 for these methods to work.")
        return

    print(f"  The effective data block size L is calculated as N - M + 1.")
    print(f"  L = {N} - {M} + 1 = {L}")
    print("-" * 50)

    # Step 3: Calculate Operations for Overlap-Add
    # The input signal (length L_x) is split into non-overlapping blocks of size L.
    # The number of operations equals the number of blocks.
    num_ops_oa = math.ceil(L_x / L)

    print("Step 3: Calculate Operations for Overlap-Add Method")
    print("  The number of operations is ceil(L_x / L).")
    print(f"  Number of operations = ceil({L_x} / {L})")
    print(f"                     = ceil({L_x / L:.3f})")
    print(f"                     = {num_ops_oa}")
    print("-" * 50)

    # Step 4: Calculate Operations for Overlap-Save
    # Each operation produces L valid output samples.
    # Total output length is L_y = L_x + M - 1.
    # The number of operations is the number of blocks to generate the full output.
    L_y = L_x + M - 1
    num_ops_os = math.ceil(L_y / L)

    print("Step 4: Calculate Operations for Overlap-Save Method")
    print(f"  First, calculate the total output length L_y = L_x + M - 1.")
    print(f"  L_y = {L_x} + {M} - 1 = {L_y}")
    print("\n  The number of operations is ceil(L_y / L).")
    print(f"  Number of operations = ceil({L_y} / {L})")
    print(f"                     = ceil({L_y / L:.3f})")
    print(f"                     = {num_ops_os}")
    print("-" * 50)

    print("Final Result:")
    print(f"For an overlap-add implementation, {num_ops_oa} (DFT+IDFT) operations are necessary.")
    print(f"For an overlap-save implementation, {num_ops_os} (DFT+IDFT) operations are necessary.")


if __name__ == '__main__':
    calculate_convolution_operations()
