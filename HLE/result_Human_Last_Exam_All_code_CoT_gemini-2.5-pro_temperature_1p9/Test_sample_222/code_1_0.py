import numpy as np

def solve_dft():
    """
    Calculates the 8-point DFT of an interleaved sequence y(n) = {x(0),h(0),...,x(3),h(3)}
    given the 4-point DFTs X(k) and H(k).
    """
    # Given 4-point DFTs
    X_k = np.array([1, 1j, -1, -1j])
    H_k = np.array([0, 1 + 1j, 1, 1 - 1j])

    print("Given 4-point DFTs:")
    print(f"X(k) = {X_k}")
    print(f"H(k) = {H_k}\n")
    print("The 8-point DFT Y(k) of the interleaved sequence is calculated using the formula:")
    print("Y(k) = X(k mod 4) + W_8^k * H(k mod 4), where W_8^k = exp(-j*2*pi*k/8)\n")
    
    # To store the results
    Y_k = np.zeros(8, dtype=np.complex128)

    print("--- Calculation for each k ---")
    for k in range(8):
        x_val = X_k[k % 4]
        h_val = H_k[k % 4]
        twiddle_factor = np.exp(-2j * np.pi * k / 8)
        y_val = x_val + twiddle_factor * h_val
        Y_k[k] = y_val

        print(f"Y({k}) = X({k % 4}) + W_8^{k} * H({k % 4})")
        print(f"   = {x_val} + ({twiddle_factor.real:.4f}{twiddle_factor.imag:+.4f}j) * {h_val}")
        print(f"   = ({y_val.real:.4f}{y_val.imag:+.4f}j)\n")

    print("\n--- Final 8-point DFT Result ---")
    print("Y(k) = [")
    for val in Y_k:
        print(f"  ({val.real:.4f}{val.imag:+.4f}j),")
    print("]")


solve_dft()

# Constructing the answer string based on the calculation.
final_answer = [
    (1.0000+0.0000j),
    (1.4142+1.0000j),
    (-1.0000-1.0000j),
    (-1.4142-1.0000j),
    (1.0000-0.0000j),
    (-1.4142+1.0000j),
    (-1.0000+1.0000j),
    (1.4142-1.0000j)
]
final_answer_str = str(final_answer).replace('j', 'j, ').replace('),', ');').replace('[', '{').replace(']', '}')

# This part is for the final required format, it won't be executed by the user.
# The code above is what the user should execute.
# print(f"<<<{final_answer_str}>>>")