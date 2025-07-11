from fractions import Fraction
import numpy as np

def solve_spin_orbit_coupling():
    """
    Calculates and prints the common eigenvalues of J^2 and J_z
    for a p-electron (l=1, s=1/2).
    """
    l = 1
    s = Fraction(1, 2)
    hbar_sq = "ħ²"
    hbar = "ħ"

    print(f"Solving for a p-electron with orbital quantum number l = {l} and spin s = {s}.\n")

    # Determine the possible values of j
    j_min = abs(l - s)
    j_max = l + s
    # Generate j values from j_min to j_max in steps of 1
    j_values = [Fraction(j) for j in np.arange(j_min, j_max + 1, 1)]

    print(f"The possible total angular momentum quantum numbers (j) range from |l-s| to l+s:")
    print(f"j = {', '.join(map(str, j_values))}\n")
    print("-" * 50)

    # Iterate through each j value to find the eigenvalues
    for j in j_values:
        print(f"For the case where j = {j}:")

        # Calculate eigenvalue of J^2
        j_plus_1 = j + 1
        j2_eigenvalue = j * j_plus_1
        print("  1. Eigenvalue of J²:")
        print(f"     Equation: j(j+1){hbar_sq}")
        print(f"     Calculation: ({j}) * ({j} + 1){hbar_sq} = ({j}) * ({j_plus_1}){hbar_sq} = {j2_eigenvalue}{hbar_sq}")

        # Determine eigenvalues of J_z
        mj_values = [Fraction(m) for m in np.arange(-j, j + 1, 1)]
        mj_str = ', '.join(map(str, mj_values))

        print("\n  2. Eigenvalues of Jz:")
        print(f"     Possible m_j values: {mj_str}")
        print(f"     Equation: m_j*{hbar}")
        print(f"     The corresponding eigenvalues are: ({mj_str})*{hbar}")
        print("-" * 50)

if __name__ == '__main__':
    solve_spin_orbit_coupling()