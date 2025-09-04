import re

def check_answer_correctness(llm_response: str) -> str:
    """
    Checks the correctness of the LLM's answer to the Maxwell's equations question.

    The function encodes the physical principles of how Maxwell's equations would change
    if magnetic monopoles existed and compares the LLM's chosen option against these principles.
    """

    # --- Define the physical principles and correct answer ---
    # 1. Gauss's Law for Magnetism (∇ ⋅ B = 0) must change to account for magnetic charge.
    #    This relates to the "divergence of the magnetic field".
    # 2. Faraday's Law (∇ × E = -∂B/∂t) must change to account for magnetic current.
    #    This relates to the "circulation of the electric field".
    # 3. The other two laws (Gauss's for E, Ampere-Maxwell) are unaffected.
    correct_changes = {"divergence of the magnetic field", "circulation of the electric field"}

    # --- Map the question's options to the physical concepts ---
    # Using simplified but descriptive keys for comparison.
    options_map = {
        "A": {"divergence of the magnetic field"},
        "B": {"divergence of the magnetic field", "curl of the magnetic field"},
        "C": {"circulation of the magnetic field", "flux of the electric field"},
        "D": {"circulation of the electric field", "divergence of the magnetic field"}
    }

    # --- Extract the LLM's chosen option ---
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Failure: Could not find a final answer in the required format '<<<X>>>'."

    chosen_option_letter = match.group(1)
    chosen_changes = options_map.get(chosen_option_letter)

    # --- Compare the chosen answer with the correct principles ---
    if chosen_changes == correct_changes:
        return "Correct"
    else:
        # --- Generate a detailed reason for the incorrectness ---
        reason_parts = []
        
        # Check for concepts that were incorrectly included
        incorrectly_added = chosen_changes - correct_changes
        if "curl of the magnetic field" in incorrectly_added or "circulation of the magnetic field" in incorrectly_added:
            reason_parts.append("it incorrectly claims the law for the circulation/curl of the magnetic field (Ampere-Maxwell Law) changes. This law is not affected by magnetic monopoles.")
        if "flux of the electric field" in incorrectly_added:
            reason_parts.append("it incorrectly claims the law for the flux of the electric field (Gauss's Law for Electricity) changes. This law is not affected by magnetic monopoles.")

        # Check for concepts that were missed
        missed_changes = correct_changes - chosen_changes
        if "divergence of the magnetic field" in missed_changes:
            reason_parts.append("it fails to identify that Gauss's Law for Magnetism (related to the divergence of the magnetic field) must change.")
        if "circulation of the electric field" in missed_changes:
            reason_parts.append("it fails to identify that Faraday's Law (related to the circulation of the electric field) must also change to account for magnetic currents.")

        # Handle the case where the answer is a correct subset (incomplete)
        if not incorrectly_added and missed_changes:
            return f"Incorrect. The answer is incomplete. While it correctly identifies one of the changing equations, it misses that the equation for the '{next(iter(missed_changes))}' also changes."

        if not reason_parts:
             return f"Incorrect. The chosen option '{chosen_option_letter}' does not correctly identify the two laws that would change."

        return f"Incorrect. The answer is wrong because {' and '.join(reason_parts)}"

# The user-provided text containing the final answer
llm_answer_text = """
Here is a step-by-step analysis to determine the correct answer.

### Step 1: Understand the core concept
The question asks how Maxwell's equations would change if magnetic monopoles (isolated North or South poles) existed. This requires understanding the physical meaning of each of the four equations in their standard form.

### Step 2: Analyze Maxwell's Equations in their standard form
1.  **Gauss's Law for Electricity (∇ ⋅ E = ρ/ε₀):** Relates the divergence of the electric field to electric charge density. It means electric field lines originate from and terminate on electric charges.
2.  **Gauss's Law for Magnetism (∇ ⋅ B = 0):** States the divergence of the magnetic field is zero. This is the crucial mathematical statement that there are no magnetic monopoles. Magnetic field lines always form closed loops.
3.  **Faraday's Law of Induction (∇ × E = -∂B/∂t):** Relates the curl (circulation) of the electric field to a changing magnetic field.
4.  **Ampère-Maxwell Law (∇ × B = μ₀(J + ε₀∂E/∂t)):** Relates the curl (circulation) of the magnetic field to electric currents and changing electric fields.

### Step 3: Determine the necessary changes for magnetic monopoles

**Careful Point 1: The direct consequence of magnetic charge.**
If a magnetic monopole (a magnetic charge) exists, it must be a source or sink for the magnetic field, just as an electric charge is for the electric field. Therefore, the divergence of the magnetic field around a monopole cannot be zero.
*   **Equation to change:** Gauss's Law for Magnetism (∇ ⋅ B = 0).
*   **New form:** It would become analogous to Gauss's Law for Electricity, looking something like ∇ ⋅ B = ρ_m, where ρ_m is the magnetic charge density.
*   **Conclusion:** The equation related to the **divergence of the magnetic field** must change.

**Careful Point 2: The consequence of magnetic current and symmetry.**
If magnetic charges can exist, they can also move, creating a "magnetic current" (J_m). The standard equations are asymmetric: an electric current (J) creates a circulating magnetic field (in the Ampère-Maxwell Law), but there is no corresponding term for a magnetic current creating an electric field. To restore this symmetry, a magnetic current term must be added.
*   **Equation to change:** Faraday's Law of Induction (∇ × E = -∂B/∂t).
*   **New form:** It would be modified to include the magnetic current as a source for a circulating electric field: ∇ × E = -∂B/∂t - J_m.
*   **Conclusion:** The equation related to the **circulation (curl) of the electric field** must also change.

**Careful Point 3: Which equations do NOT change?**
The other two equations describe how electric charges and currents create fields. The existence of magnetic monopoles does not alter these established relationships. Therefore, Gauss's Law for Electricity and the Ampère-Maxwell Law would remain unchanged.

### Step 4: Match the conclusion with the given options

Our analysis concludes that the two equations that would be different are:
1.  The one related to the **divergence of the magnetic field**.
2.  The one related to the **circulation of the electric field**.

Let's evaluate the options:
A) The one related to the divergence of the magnetic field. - *This is true but incomplete, as a second equation also changes.*
B) The ones related to the divergence and the curl of the magnetic field. - *Incorrect. The curl of the magnetic field (Ampère-Maxwell Law) does not change.*
C) The one related to the circulation of the magnetic field and the flux of the electric field. - *Incorrect. The circulation of the magnetic field (Ampère-Maxwell Law) and the flux of the electric field (Gauss's Law for Electricity) do not change.*
D) The ones related to the circulation of the electric field and the divergence of the magnetic field. - **This correctly identifies both of the modified equations.**

<<<D>>>
"""

# Run the check
result = check_answer_correctness(llm_answer_text)
print(result)