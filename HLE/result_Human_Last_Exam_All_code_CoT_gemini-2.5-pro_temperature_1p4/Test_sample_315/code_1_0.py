import textwrap

def explain_glucose_tube_experiment():
    """
    Explains the visual phenomenon of polarized white light passing through a D-glucose solution.
    """

    # --- Problem Analysis ---
    # We need to determine the visual appearance of a tube of D-glucose solution
    # when polarized white light is passed through it and viewed from the side.
    # The key physical principles are optical activity, optical rotatory dispersion,
    # and light scattering.

    # --- Explanation ---
    print("Analyzing the physics behind the observation:\n")

    # 1. Optical Activity of D-Glucose
    explanation_1 = """
    D-glucose is an optically active (chiral) molecule. This means it has the ability to rotate the plane of plane-polarized light. When polarized light enters the solution, its plane of polarization will begin to rotate as it travels.
    """
    print("Step 1: Optical Activity")
    print(textwrap.dedent(explanation_1))

    # 2. Optical Rotatory Dispersion (ORD)
    explanation_2 = """
    The angle of rotation is not the same for all colors (wavelengths) of light. This phenomenon is called Optical Rotatory Dispersion. Shorter wavelengths, like blue and violet, are rotated more per unit distance than longer wavelengths, like orange and red.
    """
    print("Step 2: Optical Rotatory Dispersion (ORD)")
    print(textwrap.dedent(explanation_2))

    # 3. Light Scattering and Observation from the Side
    explanation_3 = """
    For an observer to see the light from the side, some of it must be scattered out of the main beam by the molecules in the solution. The intensity of light scattered at 90 degrees is highly dependent on the orientation of its polarization plane relative to the observer. Maximum brightness is observed when the polarization plane is perpendicular to the line of sight.
    """
    print("Step 3: Light Scattering")
    print(textwrap.dedent(explanation_3))

    # 4. Synthesizing the Result
    explanation_4 = """
    Combining these effects: As the multi-colored polarized light travels down the tube, the polarization plane of each color rotates at a different speed. At any given point along the tube, each color has a unique polarization orientation.
    
    Therefore, at any specific point, only the color(s) with the 'correct' polarization orientation will be scattered brightly towards the observer. As you look further down the tube, the planes of all colors continue to rotate, so the color that is scattered brightly changes.
    
    This continuous change of brightly scattered color along the length of the tube, combined with the rotational nature of the phenomenon, manifests as a visible rainbow that spirals along the tube's axis.
    """
    print("Step 4: Conclusion")
    print(textwrap.dedent(explanation_4))


    # --- Final Answer ---
    answer_choice = 'C'
    answer_description = "A tube with a rainbow gradient progressing in a spiral"
    print("\n-------------------------------------------")
    print(f"The resulting appearance is therefore: '{answer_description}'")
    print(f"This corresponds to answer choice: {answer_choice}")
    print("-------------------------------------------")


if __name__ == '__main__':
    explain_glucose_tube_experiment()