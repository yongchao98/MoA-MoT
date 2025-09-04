import re

def check_symmetry_answer(llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the molecular symmetry question.

    The function uses a knowledge base of the point groups for the given molecules
    to verify the selected answer.
    """

    # --- Knowledge Base ---
    # This dictionary stores the established point group for each molecule.
    # The point group determines the molecule's symmetry properties.
    # C3h requires: 1) a C3 axis, 2) a horizontal mirror plane (σh), and 3) NO perpendicular C2 axes.
    molecule_data = {
        "triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone": {
            "point_group": "C3h",
            "reason": "This is a rigid, planar molecule. The molecular plane acts as a horizontal mirror plane (σh). It has a C3 axis perpendicular to the plane. The propeller-like fusion of the rings removes the perpendicular C2 axes that would elevate its symmetry, making its point group exactly C3h."
        },
        "benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone": {
            "point_group": "D3h",
            "reason": "This planar molecule has a C3 axis and a horizontal mirror plane (σh). However, it also possesses three C2 axes perpendicular to the C3 axis, which elevates its symmetry to the D3h point group, a higher symmetry than C3h."
        },
        "quinuclidine": {
            "point_group": "C3v",
            "reason": "This molecule has a C3 axis but lacks a horizontal mirror plane (σh). Instead, it has three vertical mirror planes (σv), placing it in the C3v point group."
        },
        "triisopropyl borate": {
            "point_group": "C3 (stable)",
            "reason": "This is a flexible molecule. While a C3h conformation is theoretically possible, its most stable, lowest-energy conformation is a twisted 'propeller' shape (C3 symmetry) to minimize steric hindrance. It is not a rigid, unambiguous example of C3h symmetry."
        }
    }

    # Mapping from the question's options to the molecule names
    options = {
        "A": "triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone",
        "B": "benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone",
        "C": "quinuclidine",
        "D": "triisopropyl borate"
    }

    # --- Answer Extraction ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<X>>> in the provided text."

    selected_option = match.group(1)
    selected_molecule_name = options.get(selected_option)

    if not selected_molecule_name:
        return f"Invalid option '{selected_option}' found."

    # --- Verification Logic ---
    # The question asks for the molecule with C3h symmetry.
    correct_molecule_name = None
    for name, data in molecule_data.items():
        if data["point_group"] == "C3h":
            correct_molecule_name = name
            break
    
    # Get the data for the selected molecule
    selected_molecule_data = molecule_data[selected_molecule_name]

    if selected_molecule_name == correct_molecule_name:
        return "Correct"
    else:
        # Construct a detailed reason for the incorrect answer
        reason = (
            f"Incorrect. The selected answer is {selected_option} ({selected_molecule_name}), "
            f"but this molecule does not have C3h symmetry.\n"
            f"Reason: {selected_molecule_data['reason']}\n"
            f"The correct molecule is '{correct_molecule_name}', which corresponds to option A. "
            f"Its point group is C3h because {molecule_data[correct_molecule_name]['reason']}"
        )
        return reason

# The final answer provided by the LLM to be checked
final_answer_text = """
Here is a step-by-step analysis to determine the correct molecule with C3h symmetry.

**1. Define the C3h Point Group**

First, it is essential to understand the symmetry elements required for a molecule to belong to the C3h point group.
*   **C3 axis**: A three-fold axis of rotation. Rotating the molecule by 120° around this axis leaves it unchanged.
*   **σh plane**: A horizontal mirror plane that is perpendicular to the C3 axis.
*   **Absence of other key elements**: A molecule is assigned to the highest symmetry point group it possesses. For a molecule to be C3h, it must lack C2 axes perpendicular to the C3 axis (which would elevate its symmetry to D3h) and it must lack vertical mirror planes (σv) if it doesn't have a σh (which would make it C3v).

**2. Analyze Each Candidate Molecule**

Let's analyze each molecule based on these requirements, using the chemical names to avoid confusion.

*   **quinuclidine**:
    *   **Structure**: A rigid, cage-like bicyclic amine.
    *   **Symmetry Analysis**: It has a C3 axis passing through the nitrogen atom and the opposite bridgehead carbon. However, it does not have a horizontal mirror plane (σh) perpendicular to this axis. Instead, it possesses three vertical mirror planes (σv) that contain the C3 axis.
    *   **Conclusion**: The presence of a C3 axis and three σv planes places quinuclidine in the **C3v** point group. It is not C3h.

*   **benzo[1,2-c:3,4-c':5,6-c'']trifuran-1,3,4,6,7,9-hexaone**:
    *   **Structure**: Also known as mellitic trianhydride, this is a planar molecule with a central benzene ring.
    *   **Symmetry Analysis**: As a planar molecule, its molecular plane acts as a horizontal mirror plane (σh). It also has a C3 axis perpendicular to this plane. However, it possesses additional symmetry: three C2 axes lie within the molecular plane, perpendicular to the C3 axis.
    *   **Conclusion**: The combination of a C3 axis, a σh plane, and three perpendicular C2 axes places this molecule in the **D3h** point group. Since D3h is a higher symmetry group than C3h, this is not the correct answer.

*   **triisopropyl borate**:
    *   **Structure**: A molecule with a central boron atom bonded to three isopropoxy groups, B(O-i-Pr)₃. This molecule is flexible due to free rotation around its single bonds.
    *   **Symmetry Analysis**: While it is theoretically possible for this molecule to adopt a conformation with C3h symmetry (with a planar B(O)₃ core acting as the σh plane), this is not its most stable state. Due to steric hindrance from the bulky isopropyl groups, the molecule's lowest-energy conformation is a twisted "propeller" shape, which belongs to the C3 point group (lacking the σh plane).
    *   **Conclusion**: Because its C3h symmetry is conformation-dependent and not inherent to a rigid, stable structure, it is a weak and ambiguous candidate.

*   **triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone**:
    *   **Structure**: A large, rigid, and planar molecule built on a triphenylene core.
    *   **Symmetry Analysis**: The molecule is planar, so the molecular plane serves as a **σh plane**. Its three-bladed shape gives it a **C3 axis** perpendicular to the plane. The crucial feature is the specific fusion pattern of the outer rings (`[1,2-c:5,6-c':9,10-c'']`). This creates a rotational, "triskelion" arrangement that breaks the symmetry of the parent triphenylene core by eliminating the three C2 axes that would be perpendicular to the C3 axis.
    *   **Conclusion**: The molecule possesses a C3 axis and a σh plane but lacks the perpendicular C2 axes. This perfectly matches the definition of the **C3h** point group, making it the best and most unambiguous answer.

**Final Decision**

After analyzing all four molecules, triphenyleno[1,2-c:5,6-c':9,10-c'']trifuran-1,3,6,8,11,13-hexaone is the only one that is a rigid molecule whose structure inherently belongs to the C3h point group.

<<<A>>>
"""

# Run the check
result = check_symmetry_answer(final_answer_text)
print(result)