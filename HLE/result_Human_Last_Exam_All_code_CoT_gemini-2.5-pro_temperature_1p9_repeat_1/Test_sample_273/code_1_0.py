def demonstrate_soliton_instability():
    """
    This function demonstrates the instability of a 3D soliton
    governed only by Heisenberg exchange and DMI terms, based on
    a scaling argument (Derrick's theorem).
    """
    print("Analyzing the stability of a hypothetical 3D soliton.")
    print("The total energy E as a function of its size scaling factor 'lambda' is:")
    print("E(lambda) = (lambda * E_exchange) + (lambda^2 * E_DMI)\n")

    # For a static (stationary) solution to exist at lambda=1, the virial
    # theorem E_exchange + 2 * E_DMI = 0 must be satisfied.
    # Since E_exchange = integral(A * (grad(m))^2) dV must be positive for any
    # localized soliton, E_DMI must be negative.
    # We assume a hypothetical solution exists that satisfies this condition.
    # The exact values don't matter, only their relationship and signs.
    E_exchange_at_1 = 2.0
    E_DMI_at_1 = -1.0

    print("Let's assume a hypothetical stationary soliton exists at lambda=1 with:")
    print(f"  - Exchange Energy at lambda=1 (E_exchange) = {E_exchange_at_1}")
    print(f"  - DMI Energy at lambda=1 (E_DMI)        = {E_DMI_at_1}")
    
    # Check the stationary point condition
    print("\nThis choice satisfies the condition for a stationary point, as:")
    equation_sum = E_exchange_at_1 + 2 * E_DMI_at_1
    # We output each number in the equation.
    print(f"  dE/d(lambda)|_1 = E_exchange + 2 * E_DMI = {E_exchange_at_1} + 2 * ({E_DMI_at_1}) = {equation_sum}\n")

    # Now, let's check for stability by examining the second derivative.
    # Stability requires the second derivative of E(lambda) to be positive at lambda=1.
    # d^2(E)/d(lambda^2) = 2 * E_DMI
    second_derivative = 2 * E_DMI_at_1
    print("For the soliton to be stable, the energy at lambda=1 must be a minimum.")
    print("This means the second derivative of the energy must be positive.")
    # We output each number in the equation.
    print(f"  d^2(E)/d(lambda^2)|_1 = 2 * E_DMI = 2 * ({E_DMI_at_1}) = {second_derivative}")

    if second_derivative > 0:
        print("--> The second derivative is positive, indicating stability (a MINIMUM).")
    else:
        print("--> The second derivative is negative, indicating the solution is UNSTABLE (a MAXIMUM).\n")

    # Let's explicitly calculate the total energy for different sizes (lambda) to verify.
    print("We can verify this by calculating the total energy for different soliton sizes:")
    
    lambdas_to_check = [0.8, 0.9, 1.0, 1.1, 1.2]
    
    print("-" * 62)
    print(" Lambda | Eq: (l*E_ex) + (l^2*E_DMI)           | Total Energy")
    print("-" * 62)

    for l in lambdas_to_check:
        # E(lambda) = lambda * E_ex(1) + lambda^2 * E_DMI(1)
        exchange_term = l * E_exchange_at_1
        dmi_term = (l**2) * E_DMI_at_1
        total_energy = exchange_term + dmi_term
        
        # Outputting each number in the final equation as requested
        equation_str = f"({l:.1f}*{E_exchange_at_1}) + ({l:.1f}²*{E_DMI_at_1})"
        
        status = ""
        if l == 1.0:
            status = "<-- Unstable Stationary Point"

        print(f" {l:^6.1f} | {equation_str:<32} | {total_energy:^12.4f} {status}")

    print("-" * 62)
    print("\nConclusion: The energy is maximum at the stationary point (lambda=1).")
    print("Any small change in size (shrinking for lambda<1 or expanding for lambda>1)")
    print("lowers the total energy. The soliton is therefore fundamentally unstable.")


demonstrate_soliton_instability()