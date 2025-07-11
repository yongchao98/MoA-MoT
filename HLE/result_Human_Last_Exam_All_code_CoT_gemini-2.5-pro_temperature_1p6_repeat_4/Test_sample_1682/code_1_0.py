def solve_quiz():
    """
    Analyzes historical facts about Chinese wedding attire to find the incorrect statement.
    """
    incorrect_statement = "C"
    explanation = [
        "The incorrect statement is C.",
        "Explanation: While it was customary for brides in the Tang and Song dynasties to use a fan to cover their faces (a practice known as 以扇遮面, yǐ shàn zhē miàn), the type of fan specified is not entirely accurate for the entire period.",
        "Folding fans (摺扇, zhéshàn) were introduced to China and only became widely popular during the Song Dynasty.",
        "In the earlier Tang Dynasty, brides would have used a round fan (團扇, tuánshàn).",
        "Therefore, the statement that it was customary to use *folding fans* during the Tang Dynasty is incorrect."
    ]

    for line in explanation:
        print(line)

    print(f"<<<{incorrect_statement}>>>")

solve_quiz()