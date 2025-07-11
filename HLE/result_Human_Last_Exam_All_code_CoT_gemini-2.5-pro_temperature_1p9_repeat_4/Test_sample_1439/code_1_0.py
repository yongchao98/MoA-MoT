import sys

def solve_critical_exponent_order():
    """
    This script explains the order in the coupling constant 'u' at which the
    critical exponent ν acquires its first non-vanishing contribution in the
    perturbative epsilon-expansion of φ⁴ theory.
    """
    
    print("Step-by-step analysis of the critical exponent ν in φ⁴ theory:")
    print("-" * 75)

    # Step 1: Mean-Field Theory Baseline
    nu_0 = 1 / 2
    print(f"1. In mean-field theory (corresponding to a non-interacting system with coupling u = 0), "
          f"the critical exponent ν has the value ν₀ = {nu_0}.")

    # Step 2: Introducing Interactions
    print("\n2. The Renormalization Group (RG) method allows for systematic calculation of corrections to this value "
          "due to interactions, which are described by a coupling constant 'u'.")

    # Step 3: Perturbative Expansion
    print("\n3. Within the ε-expansion framework, physical quantities are expressed as a power series "
          "in the coupling constant 'u' evaluated at the non-trivial (Wilson-Fisher) fixed point, u*.")
    
    # Step 4: Relation between ν and anomalous dimension
    print("\n4. The exponent ν is not calculated directly, but is derived from other quantities. Specifically, it is related "
          "to an 'anomalous dimension', γ_m(u), associated with the mass term in the theory.")
    print("   The relationship is: 1/ν = 2 + γ_m(u).")

    # Step 5: Perturbative nature of the anomalous dimension
    print("\n5. The anomalous dimension γ_m(u) is calculated using Feynman diagrams as a power series in u. "
          "The crucial point is that the lowest-order (one-loop) diagrammatic contribution to γ_m(u) "
          "is non-zero and is proportional to the coupling 'u' itself.")
    print("   This means the expansion for γ_m(u) starts at the first power of u:")
    print("   γ_m(u) = A * u¹ + B * u² + ... (where A and B are constants)")

    # Step 6: Deriving the order for ν
    print("\n6. We can now see how the correction to ν appears. By substituting the series for γ_m(u) back into the equation for ν:")
    print("   ν = 1 / (2 + A*u + ...)")
    print("   Expanding this for small u (using the Taylor series (1+x)⁻¹ ≈ 1-x) gives:")

    nu_zero_val = 0.5
    order_u = 1
    
    print("\n   ν(u) ≈ (1/2) * (1 - (A/2)*u)")
    print(f"   ν(u) = {nu_zero_val} + c₁*u^{order_u} + O(u²)")
    print("   where c₁ is a non-zero constant.")
    print("-" * 75)
    
    # Final Conclusion
    print("\nConclusion:")
    print("The equation shows that the first term that modifies ν from its mean-field value of 1/2 is the term linear in 'u'.")
    print("Therefore, the initial non-vanishing contribution to the critical exponent ν appears at the first order in the coupling constant u.")

if __name__ == '__main__':
    solve_critical_exponent_order()