import re
import sympy

def check_physics_answer(llm_response_text):
    """
    Checks the correctness of the LLM's answer to the electrostatics problem.

    The function verifies the answer based on the principle of electrostatic shielding,
    which dictates that the external electric field depends only on the total charge
    enclosed and the distance from the center of the spherical conductor (L), not on
    the internal geometry (l, s, theta).
    """
    # Define symbolic variables to represent the physical quantities
    q, L, l, s, theta, epsilon_o = sympy.symbols('q L l s theta epsilon_o')
    k = 1 / (4 * sympy.pi * epsilon_o)

    # Define the options as presented in the LLM's final analysis
    # This mapping is consistent across most of the provided candidate answers
    options = {
        'A': k * q / (l - s * sympy.cos(theta))**2,
        'B': k * q / L**2,
        'C': k * q / l**2,
        'D': k * q / (l + s * sympy.cos(theta))**2
    }

    # Principle 1: The correct formula must be independent of internal geometry (l, s, theta)
    # due to electrostatic shielding.
    # Principle 2: The field should follow the inverse square law with respect to the distance
    # from the center of the conductor (L), according to the Shell Theorem.
    correct_formula = k * q / L**2
    
    correct_option_letter = None
    for letter, formula in options.items():
        if formula == correct_formula:
            correct_option_letter = letter
            break
            
    if correct_option_letter is None:
        # This case should not happen if the options are as described
        return "Error in problem analysis: The theoretically correct formula is not listed among the options."

    # Extract the final answer choice from the LLM's response
    match = re.search(r'<<<([A-D])>>>\s*$', llm_response_text, re.MULTILINE)
    if not match:
        return "Could not find the final answer in the required format '<<<X>>>' at the end of the response."
        
    llm_choice = match.group(1)
    
    # Compare the LLM's choice with the physically correct option
    if llm_choice == correct_option_letter:
        return "Correct"
    else:
        chosen_formula = options.get(llm_choice)
        
        reason = (f"Incorrect. The final answer was given as option {llm_choice}, "
                  f"which corresponds to the formula E = {chosen_formula}.\n\n")
        
        # Analyze why the chosen formula is wrong
        free_symbols = chosen_formula.free_symbols
        if l in free_symbols or s in free_symbols or theta in free_symbols:
            reason += ("This is incorrect because the principle of electrostatic shielding dictates that the "
                       "electric field outside the conductor cannot depend on the internal geometry "
                       "of the cavity (variables l, s, or theta).")
        elif L not in free_symbols:
             reason += ("This is incorrect because the field must depend on the distance L from the "
                        "center of the conductor.")
        else:
             reason += "This formula has an incorrect form."

        reason += (f"\n\nThe correct answer is option {correct_option_letter}, E = {correct_formula}, "
                   "because the charge +q on the outer surface distributes uniformly, and the field "
                   "outside is equivalent to that of a point charge at the conductor's center.")
        
        return reason

# The user-provided text containing the question and the final answer from the LLM
llm_full_response = """
Here is a step-by-step analysis of the problem to determine the correct answer.

**1. Analyze the Initial State and Charge Induction:**
*   The system starts with an uncharged spherical conductor. This means its net charge is zero.
*   A positive point charge `+q` is placed inside the cavity.
*   Due to electrostatic induction, this positive charge attracts the free electrons within the conductor. These electrons move and accumulate on the surface of the cavity.
*   To determine the amount of induced charge on the inner surface, we use Gauss's Law. Imagine a Gaussian surface entirely within the conducting material, enclosing the cavity. In electrostatic equilibrium, the electric field inside the conductor's material is zero (`E=0`).
*   Since the flux (`∮ E ⋅ dA`) through this Gaussian surface is zero, the net charge enclosed (`Q_enclosed`) must also be zero.
*   The enclosed charge is the sum of the point charge `+q` and the charge induced on the inner surface, `q_inner`. So, `+q + q_inner = 0`, which means `q_inner = -q`.
*   A charge of `-q` is induced on the inner surface of the cavity.

**2. Conservation of Charge and the Outer Surface:**
*   The conductor was initially neutral. Since a charge of `-q` has moved to the inner surface, a charge of `+q` must be left on the outer surface of the conductor to maintain overall charge neutrality. Let's call this `q_outer = +q`.

**3. The Principle of Electrostatic Shielding:**
*   This is the most crucial concept for solving the problem. A conducting shell acts as a "Faraday cage". It shields the region outside the conductor from the electric fields of any charges located inside it.
*   The electric field at any point P outside the conductor is the superposition of the fields from all charges: the original `+q`, the induced `q_inner = -q`, and the induced `q_outer = +q`.
*   The shielding principle tells us that the combined electric field from the charge inside the cavity (`+q`) and the charge on the inner surface (`-q`) is exactly zero everywhere outside the cavity.
*   Therefore, the electric field at point P is determined *solely* by the charge `q_outer = +q` on the outer surface of the conductor.

**4. Distribution of Charge on the Outer Surface:**
*   The charge `+q` resides on the outer surface of the conductor.
*   Because the outer surface is a sphere and it is shielded from the asymmetric arrangement of charges inside, this charge `+q` will distribute itself *uniformly* over the spherical surface. This is the lowest energy configuration and ensures the outer surface is an equipotential.

**5. Calculating the Electric Field at Point P:**
*   We need to find the electric field at point P due to a uniform spherical shell of charge `+q` and radius `R`.
*   According to the Shell Theorem (a consequence of Gauss's Law), the electric field at any point outside a uniformly charged spherical shell is identical to the electric field that would be produced if all the charge were concentrated at the center of the sphere.
*   In this case, we can treat the `+q` on the outer surface as a point charge located at the center of the spherical conductor.
*   The distance from the center of the conductor to point P is given as `L`.
*   Using Coulomb's Law for a point charge, the magnitude of the electric field is:
    $E = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{\\text{charge}}{\\text{distance}^2} = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{q}{L^2}$

**Careful Points and Conclusion:**
*   The information about the cavity's radius (`r`), its position (`s`), the distance from the cavity's center to P (`l`), and the angle (`θ`) are all distractors. They describe the internal geometry, which is irrelevant to the external electric field due to electrostatic shielding.
*   The final electric field depends only on the total charge placed inside the cavity (`q`) and the distance from the center of the *outer conductor* to the point of measurement (`L`).
*   Comparing our result, $E = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{q}{L^2}$, with the given options:
    A) $E = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{q}{(l-s \\cos\\theta)^2}$
    B) $E = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{q}{L^2}$
    C) $E = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{q}{l^2}$
    D) $E = \\dfrac{1}{4 \\pi \\epsilon_o} \\dfrac{q}{(l+s \\cos\\theta)^2}$
*   The correct option is B.

<<<B>>>
"""

# Run the check
result = check_physics_answer(llm_full_response)
print(result)