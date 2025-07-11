# This script provides an explanation for the user's question about
# spontaneous parametric down-conversion (SPDC) in boron nanosheets.
# The explanation is stored in a string and then printed to the console.

def explain_spdc_in_borophene():
    """
    This function prints a detailed answer to the question of whether
    free-standing boron nanosheets exhibit SPDC.
    """
    # The answer hinges on a key principle of nonlinear optics: material symmetry.
    # SPDC is a second-order nonlinear process.
    principle = "Spontaneous Parametric Down-Conversion (SPDC) is a second-order nonlinear optical process governed by the material's χ⁽²⁾ susceptibility."

    # Such processes have a strict requirement for the material's crystal structure.
    requirement = "A fundamental requirement for a material to exhibit a non-zero bulk χ⁽²⁾ response is that its crystal structure must lack a center of inversion (i.e., it must be non-centrosymmetric)."

    # Now, we apply this to the specific material: boron nanosheets.
    material_properties = "The most common and stable phases of free-standing boron nanosheets (borophene), such as the β₁₂ and χ₃ phases, are known to be centrosymmetric. They possess a center of inversion symmetry."

    # The conclusion follows directly from applying the principle to the material.
    conclusion = "Because their crystal structure has inversion symmetry, their bulk second-order susceptibility χ⁽²⁾ is zero. Therefore, second-order processes like SPDC are forbidden by this symmetry."

    # Final summary statement.
    final_answer = "In conclusion, free-standing boron nanosheets would not be expected to exhibit spontaneous parametric down-conversion."

    # Printing the full explanation. There is no equation with numbers in this problem.
    print(principle)
    print(requirement)
    print(material_properties)
    print(conclusion)
    print(final_answer)

# Execute the function to provide the answer.
explain_spdc_in_borophene()