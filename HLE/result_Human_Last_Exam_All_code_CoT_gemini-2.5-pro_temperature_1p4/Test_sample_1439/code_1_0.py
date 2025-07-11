def solve_critical_exponent_order():
    """
    This script explains and determines the order in the coupling constant u
    at which the critical exponent ν first receives a non-vanishing contribution
    in the RG analysis of φ⁴ theory.
    """

    # Step 1: Explain the theoretical background.
    # The critical exponent ν is related to the coupling u via an exact formula from the Renormalization Group.
    # We will represent this relationship and then perform a perturbative expansion.
    explanation = """
Step 1: Theoretical Background

In the renormalization group (RG) analysis of φ⁴ theory, the critical exponent ν, which governs the divergence of the correlation length, can be expressed as a function of the interaction coupling constant `u`. The exact relation is given by:

ν(u) = 1 / (2 - γ_φ²(u))

where γ_φ²(u) is the anomalous dimension of the φ² operator. In the mean-field approximation (which corresponds to the non-interacting case where `u=0`), γ_φ²(0) = 0, leading to the classical value ν = 1/2.

Step 2: Perturbative Expansion

To find corrections to the mean-field value, γ_φ²(u) is calculated as a power series in `u` by evaluating Feynman diagrams. The lowest-order (one-loop) diagram contributing to the renormalization of the φ² operator yields a term that is linear in `u`. Thus, the series expansion for γ_φ²(u) starts as:

γ_φ²(u) = c₁*u + c₂*u² + O(u³)

where c₁ is a non-zero constant.

Step 3: Finding the First Contribution to ν

We can substitute this series into the expression for ν(u) and expand it for small `u`:

ν(u) = 1 / (2 - (c₁*u + O(u²)))
     = (1/2) * [1 / (1 - (c₁/2)*u + O(u²))]

Using the geometric series expansion 1/(1-x) ≈ 1 + x for small x, we get:

ν(u) ≈ (1/2) * [1 + (c₁/2)*u]
     = 1/2 + (c₁/4)*u + O(u²)

Step 4: Conclusion

The expression `ν(u) = 1/2 + (c₁/4)*u + O(u²)` shows the value of ν as a series in the coupling constant `u`.
- The first term, `1/2`, is the mean-field value (order u⁰).
- The second term, `(c₁/4)*u`, is the first correction to the mean-field value.

This correction term is of order u¹. Therefore, the initial non-vanishing contribution to the critical exponent ν appears at the first order in the coupling constant `u`.
"""
    print(explanation)

    # The final answer is the order of the first non-vanishing contribution.
    order = 1
    
    # Per the instruction "output each number in the final equation",
    # we can represent the final answer as the equation for the order.
    print("The final equation for the order is:")
    print(f"order = {order}")


solve_critical_exponent_order()