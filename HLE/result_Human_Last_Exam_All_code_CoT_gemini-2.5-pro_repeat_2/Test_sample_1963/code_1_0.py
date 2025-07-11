import numpy as np

def solve_bootstrap():
    """
    Performs a grid search for the minimal E and <x^2> for the potential V(x) = x^2 + x^4
    using the quantum bootstrap method with K=7.
    """
    # Set a small tolerance for checking non-negativity of eigenvalues
    tolerance = -1e-9

    # Search range for Energy (E) and <x^2>
    # Literature suggests E_0 is around 1.392, so we search around this value.
    e_range = np.arange(1.3, 1.5, 0.001)
    x2_range = np.arange(0.3, 0.5, 0.001)

    for E_val in e_range:
        for x2_val in x2_range:
            # Y_k stores <x^(2k)>. We need moments up to Y_7 = <x^14>
            # because M_ij = <x^(i+j)> and for K=7, max(i+j)=14.
            Y = np.zeros(8)
            
            # Step 1: Initialize the first moments
            Y[0] = 1.0  # <x^0> = 1
            Y[1] = x2_val  # <x^2> is a search parameter

            # Step 2: Calculate higher moments using the recursion relations
            # Relation from t=1: E = 2*<x^2> + 3*<x^4>  => E = 2*Y_1 + 3*Y_2
            Y[2] = (E_val - 2 * Y[1]) / 3.0

            # Relation from t=3: 6 + 12*E*Y_1 - 16*Y_2 - 20*Y_3 = 0
            Y[3] = (6.0 + 12.0 * E_val * Y[1] - 16.0 * Y[2]) / 20.0

            # General recursion for Y_{m+2} derived from the problem statement for m>=2
            # (8m+12)Y_{m+2} = (2m+1)(2m)(2m-1)Y_{m-1} + 4(2m+1)E Y_m - (8m+8)Y_{m+1}
            for m in range(2, 6):
                # We need to calculate Y[4], Y[5], Y[6], Y[7]
                # m=2 -> Y[4]
                # m=3 -> Y[5]
                # m=4 -> Y[6]
                # m=5 -> Y[7]
                numerator = ((2*m+1)*(2*m)*(2*m-1) * Y[m-1] +
                             4*(2*m+1)*E_val * Y[m] -
                             (8*m+8) * Y[m+1])
                denominator = (8*m+12)
                if abs(denominator) < 1e-12:
                    # Should not happen in this range
                    continue
                Y[m+2] = numerator / denominator

            # Step 3: Form the positive semidefinite matrices
            # Matrix for even powers in the operator O
            M_even = np.array([
                [Y[0], Y[1], Y[2], Y[3]],
                [Y[1], Y[2], Y[3], Y[4]],
                [Y[2], Y[3], Y[4], Y[5]],
                [Y[3], Y[4], Y[5], Y[6]]
            ])

            # Matrix for odd powers in the operator O
            M_odd = np.array([
                [Y[1], Y[2], Y[3], Y[4]],
                [Y[2], Y[3], Y[4], Y[5]],
                [Y[3], Y[4], Y[5], Y[6]],
                [Y[4], Y[5], Y[6], Y[7]]
            ])

            # Step 4: Check the positive semidefinite condition
            try:
                eig_even = np.linalg.eigvalsh(M_even)
                eig_odd = np.linalg.eigvalsh(M_odd)
            except np.linalg.LinAlgError:
                continue # Matrix was not well-formed, skip

            if np.all(eig_even >= tolerance) and np.all(eig_odd >= tolerance):
                # Found the minimal E and corresponding <x^2>
                print(f"Minimal value of E: {E_val:.3f}")
                print(f"Minimal value of <x^2>: {x2_val:.3f}")
                
                # Outputting the numbers for the final matrix equations as requested
                print("\nThe corresponding moments Y_k = <x^(2k)> for k=0 to 7 are:")
                moment_strs = [f"{y:.4f}" for y in Y]
                print(moment_strs)
                
                # The problem asks for minimal E and <x^2> values
                # We format the final answer as requested
                final_answer = f"E={E_val:.3f}, <x^2>={x2_val:.3f}"
                print(f"\n<<<{final_answer}>>>")
                return

    print("No solution found in the specified range.")

if __name__ == '__main__':
    solve_bootstrap()