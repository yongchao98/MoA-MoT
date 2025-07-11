import textwrap

def evaluate_teacher_judgment_statements():
    """
    Analyzes statements about teacher judgment accuracy based on empirical research
    to identify the incorrect one.
    """

    statements = [
        {
            "option": "A",
            "statement": "Experienced teachers often overestimate students' academic achievement.",
            "is_correct": False,
            "justification": "This statement is not consistently correct. While some meta-analyses show a slight average tendency to overestimate, teacher judgment bias is complex. Teachers may also underestimate certain groups of students (e.g., high-achievers or those from lower socio-economic backgrounds). Therefore, stating they 'often' overestimate is an oversimplification and not a robustly supported finding compared to the other statements."
        },
        {
            "option": "B",
            "statement": "Experienced teacher judgment accuracy correlates with students' achievement around r = .6",
            "is_correct": True,
            "justification": "This is a well-established finding. For example, a major meta-analysis by Südkamp, Kaiser, & Möller (2012) found an average correlation of r = .63 between teacher judgments and student achievement scores. Thus, r = .6 is a commonly cited and accurate approximation."
        },
        {
            "option": "C",
            "statement": "Experienced teachers have more difficulty accurately judging low-performing students.",
            "is_correct": True,
            "justification": "This is supported by research. Judgment accuracy tends to be lower for students at the extremes of the performance distribution (both low- and high-achievers) compared to students in the middle. Difficulties in judging low-performers is a consistent finding."
        },
        {
            "option": "D",
            "statement": "Experienced teachers are typically more accurate in assessing student performance than student teachers.",
            "is_correct": True,
            "justification": "Research on expertise supports this statement. Experienced teachers generally demonstrate higher, though not overwhelmingly superior, judgment accuracy compared to pre-service or novice teachers."
        },
        {
            "option": "E",
            "statement": "Experienced teachers’ judgement accuracy on students achievement is higher than on students motivation.",
            "is_correct": True,
            "justification": "This is a consistent finding. Teachers are generally more accurate when judging cognitive characteristics like achievement (which are often directly assessed) than when judging non-cognitive or affective characteristics like motivation, which are more latent and harder to observe."
        }
    ]

    incorrect_statement = None
    for item in statements:
        if not item['is_correct']:
            incorrect_statement = item
            break

    if incorrect_statement:
        print("The incorrect statement is:\n")
        print(f"Option {incorrect_statement['option']}: {incorrect_statement['statement']}\n")
        print("Justification:")
        print(textwrap.fill(incorrect_statement['justification'], 80))
    else:
        print("Could not identify an incorrect statement based on the provided research data.")

evaluate_teacher_judgment_statements()
<<<A>>>