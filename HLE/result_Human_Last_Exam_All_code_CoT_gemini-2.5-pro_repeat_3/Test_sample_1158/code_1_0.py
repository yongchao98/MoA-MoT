import textwrap

def solve_oled_problem():
    """
    This function analyzes the properties of air-stable organic radicals in OLEDs
    to determine their main disadvantage and select the correct answer from a list.
    """

    explanation = """
    The question asks for the primary disadvantage of air-stable organic radicals in OLEDs. Let's analyze the options:

    1.  **Premise:** Radicals in OLEDs emit from doublet states, which allows them to bypass the spin statistics limitations of traditional fluorescent (singlet) and phosphorescent (triplet) materials, theoretically enabling 100% internal quantum efficiency (IQE).

    2.  **Analysis of Choices:**
        *   **A. not stable to oxygen because it is highly reactive:** This is incorrect. The question explicitly specifies "air-stable" radicals, meaning they are designed to be stable in the presence of oxygen.
        *   **B. wide FWDH because it has multiple emissions:** While color purity (related to FWHM - Full Width at Half Maximum) is important, it is not the most fundamental and universal disadvantage for this entire class of materials.
        *   **C. low luminance because radical can quench with each other:** This is a correct observation. High concentrations of radicals in the emissive layer can lead to self-quenching (or annihilation), where excited radicals interact and decay non-radiatively. This reduces light output (luminance).
        *   **D. low EQE because excitons can be quenched by the radicals:** This is the most fundamental explanation. The "exciton" in a radical-based OLED is the excited radical itself. The quenching process mentioned in (C) is a non-radiative decay pathway. This directly lowers the efficiency of converting excitons into photons, which is measured by the External Quantum Efficiency (EQE). Low luminance is a consequence of low EQE. Therefore, this option describes the root cause more accurately than (C).
        *   **E. low EQE because delocalization of the radical will have multiple emissions:** This is incorrect. Radical delocalization is a key feature for achieving stability and does not inherently cause multiple emissions that would lead to a low EQE.

    3.  **Conclusion:** The primary disadvantage is the intermolecular quenching process that provides a strong non-radiative decay channel, fundamentally limiting the device's efficiency. Option D provides the best description of this issue.
    """

    print(textwrap.dedent(explanation).strip())

    # The final answer is D
    final_answer = 'D'
    print(f"\nFinal Answer: <<<D>>>")

solve_oled_problem()