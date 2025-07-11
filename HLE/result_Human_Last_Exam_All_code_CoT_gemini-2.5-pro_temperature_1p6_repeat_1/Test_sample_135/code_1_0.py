def find_incorrect_statement():
    """
    This function models a knowledge base of empirical research on teacher judgment
    accuracy and identifies which of the provided statements is not correct.
    """
    research_findings = [
        {
            "option": "A",
            "statement": "Experienced teachers often overestimate students' academic achievement.",
            "is_correct": True,
            "reasoning": "This is a well-documented leniency or positivity bias in teacher judgments."
        },
        {
            "option": "B",
            "statement": "Experienced teacher judgment accuracy correlates with students' achievement around r = .6",
            "is_correct": True,
            "reasoning": "Meta-analyses (e.g., Hattie, 2009; Südkamp et al., 2012) report average correlations around r = .63, making this a correct statement."
        },
        {
            "option": "C",
            "statement": "Experienced teachers have more difficulty accurately judging low-performing students.",
            "is_correct": False,
            "reasoning": "Research often shows teachers are more accurate at judging students at the performance extremes (high and low) than students in the middle. Therefore, this statement is generally incorrect."
        },
        {
            "option": "D",
            "statement": "Experienced teachers are typically more accurate in assessing student performance than student teachers.",
            "is_correct": True,
            "reasoning": "Expertise and experience are positively correlated with judgment accuracy."
        },
        {
            "option": "E",
            "statement": "Experienced teachers’ judgement accuracy on students achievement is higher than on students motivation.",
            "is_correct": True,
            "reasoning": "Teachers are significantly more accurate when judging observable cognitive traits (achievement) than less observable affective traits (motivation)."
        }
    ]

    incorrect_statement = None
    for finding in research_findings:
        if not finding["is_correct"]:
            incorrect_statement = finding
            break

    if incorrect_statement:
        print("The incorrect statement is:")
        print(f"{incorrect_statement['statement']}")
    else:
        print("All statements were found to be correct according to the knowledge base.")

find_incorrect_statement()