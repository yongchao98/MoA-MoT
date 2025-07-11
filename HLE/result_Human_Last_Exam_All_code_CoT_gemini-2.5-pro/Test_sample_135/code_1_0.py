def analyze_teacher_judgment_statements():
    """
    Analyzes statements about teacher judgment accuracy to find the incorrect one.
    """
    statements = {
        'A': {
            'text': "Experienced teachers often overestimate students' academic achievement.",
            'is_correct': True,
            'reason': "Meta-analyses show a general tendency for teachers to be overly optimistic about student performance, a finding sometimes referred to as a 'positivity bias'."
        },
        'B': {
            'text': "Experienced teacher judgment accuracy correlates with students' achievement around r = .6",
            'is_correct': True,
            'reason': "This is a well-established finding. Major meta-analyses (e.g., Südkamp, Kaiser, & Möller, 2012) report an average correlation between teacher judgments and standardized test scores of r = .63, making r = .6 a correct and representative value."
        },
        'C': {
            'text': "Experienced teachers have more difficulty accurately judging low-performing students.",
            'is_correct': False,
            'reason': "This statement is incorrect. Research indicates that teachers are generally more accurate at identifying students at the extremes of the performance spectrum (both low- and high-achievers) than they are at differentiating among the large group of average-performing students."
        },
        'D': {
            'text': "Experienced teachers are typically more accurate in assessing student performance than student teachers.",
            'is_correct': True,
            'reason': "Experience generally correlates with judgment accuracy. Experienced teachers have more developed internal standards and reference points, leading to more accurate assessments than those of novices or student teachers."
        },
        'E': {
            'text': "Experienced teachers’ judgement accuracy on students achievement is higher than on students motivation.",
            'is_correct': True,
            'reason': "Judgments are more accurate for concrete, directly observable constructs like achievement (verified by tests/assignments) than for abstract, internal constructs like motivation or personality."
        }
    }

    incorrect_statement_key = None
    for key, value in statements.items():
        if not value['is_correct']:
            incorrect_statement_key = key
            break

    print("Analysis of statements on teacher judgment accuracy:")
    for key, value in statements.items():
        status = "CORRECT" if value['is_correct'] else "INCORRECT"
        print(f"\nStatement {key}: \"{value['text']}\"")
        print(f"  - Evaluation: This statement is {status} based on empirical research.")
        print(f"  - Rationale: {value['reason']}")

    if incorrect_statement_key:
        print("\n" + "="*50)
        print("CONCLUSION:")
        print(f"The statement that is NOT correct is option '{incorrect_statement_key}'.")
        print(f"It incorrectly claims that teachers have more difficulty judging low-performing students.")
        print("="*50)
    else:
        print("\nCould not identify an incorrect statement.")

if __name__ == '__main__':
    analyze_teacher_judgment_statements()
<<<C>>>