def solve_partition_function():
    """
    Derives the partition function Z for a system with Hamiltonian H = -μN
    in the grand canonical ensemble using a path integral framework.
    """
    print("This script derives the partition function Z for the given system.")
    print("="*60)

    print("\nStep 1: State the Grand Canonical Partition Function Z")
    print("The general definition of the partition function in the grand canonical ensemble is:")
    print("Z = Tr(exp(-β * (H - μ*N)))")
    print("Here, β is the inverse temperature (1/kT), H is the Hamiltonian, μ is the chemical potential, and N is the particle number operator.")

    print("\nStep 2: Substitute the System's Hamiltonian")
    print("The problem provides the Hamiltonian H = -μ*N.")
    print("Substituting this into the general formula:")
    print("Z = Tr(exp(-β * ((-μ*N) - μ*N)))")
    print("Z = Tr(exp(-β * (-2*μ*N)))")
    print("Z = Tr(exp(2*β*μ*N))")

    print("\nStep 3: Interpret the Trace as a Path Integral (Sum over Histories)")
    print("A path integral sums over all possible histories. The operator in the exponent, 2*β*μ*N, is diagonal in the particle number basis |n>.")
    print("This means the particle number 'n' is a conserved quantity, so the only valid 'paths' are those with a constant particle number.")
    print("Therefore, the trace (sum over states) becomes a simple sum over all possible particle numbers n = 0, 1, 2, ...")

    print("\nStep 4: Determine the Contribution of Each Path")
    print("The contribution of each path 'n' is exp(-S_n), where S_n is the Euclidean action.")
    print("The action is the integral of the energy eigenvalue over the imaginary time interval [0, β].")
    print("The energy operator for our expression is K = -2*μ*N, with eigenvalues E_n = -2*μ*n.")
    print("S_n = Integral from 0 to β of (E_n) dτ = Integral(-2*μ*n) dτ = -2*β*μ*n")
    print("The contribution for path 'n' is exp(-S_n) = exp(-(-2*β*μ*n)) = exp(2*β*μ*n).")

    print("\nStep 5: Sum Over All Paths to Get Z")
    print("Z is the sum of contributions from all possible particle numbers:")
    print("Z = Σ_{n=0 to ∞} exp(2*β*μ*n)")

    print("\nStep 6: Evaluate the Geometric Series")
    print("This is a geometric series of the form 1 + x + x^2 + ... with x = exp(2*β*μ).")
    print("The series converges to 1 / (1 - x) if |x| < 1, which implies μ < 0.")
    print("Assuming convergence, the partition function is:")
    print("Z = 1 / (1 - exp(2*β*μ))")

    print("\n" + "="*60)
    print("FINAL RESULT")
    print("="*60)
    print("The final equation for the partition function is: Z = 1 / (1 - exp(2*β*μ))")
    print("\nBreaking down the numbers and symbols in the final equation:")
    print("Numerator: 1")
    print("First term in denominator: 1")
    print("Operator in denominator: -")
    print("Function in denominator: exp( ... )")
    print("Argument of exp: 2 * β * μ")

if __name__ == "__main__":
    solve_partition_function()