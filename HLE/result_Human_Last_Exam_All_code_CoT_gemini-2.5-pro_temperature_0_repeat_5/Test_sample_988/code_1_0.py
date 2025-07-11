def solve_antioxidant_question():
    """
    This function identifies and explains the initial antioxidant response in
    Microcystis aeruginosa to high temperature stress.
    """
    question = "Which antioxidants are initially activated in Microcystis aeruginosa CAAT 2005-3 to counteract oxidative stress in response to high temperature exposure (29ºC)?"

    answer_choices = {
        'A': "Liposoluble antioxidants",
        'B': "Hydrosoluble antioxidants",
        'C': "Enzymatic antioxidants",
        'D': "Photosynthetic pigments",
        'E': "UV-protective compounds"
    }

    correct_answer_key = 'C'

    explanation = """
High temperature exposure induces the production of Reactive Oxygen Species (ROS) in photosynthetic organisms like Microcystis aeruginosa, leading to oxidative stress. The cell's first and most rapid line of defense is the enzymatic antioxidant system.

1.  Superoxide Dismutase (SOD), a key enzyme, is immediately activated to dismutate the highly reactive superoxide radical (O₂⁻) into hydrogen peroxide (H₂O₂).
2.  Following this, other enzymes like Catalase (CAT) and Ascorbate Peroxidase (APX) are activated to detoxify the H₂O₂.

This enzymatic response is the most immediate, preceding the slower synthesis or mobilization of non-enzymatic antioxidants (liposoluble or hydrosoluble). Therefore, enzymatic antioxidants are the ones *initially* activated.
"""

    print(explanation)
    print("Question: " + question)
    print("\nAnswer Choices:")
    for key, value in answer_choices.items():
        print(f"{key}. {value}")

    print("\n----------------------------------")
    print("Final Answer:")
    print(f"{correct_answer_key}. {answer_choices[correct_answer_key]}")
    print("----------------------------------")

solve_antioxidant_question()