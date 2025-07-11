import textwrap

def analyze_teacher_judgment_statements():
    """
    Analyzes statements about empirical research on teacher judgment accuracy
    to identify the incorrect one.
    """
    
    analysis = {
        'A': {
            'statement': "Experienced teachers often overestimate students' academic achievement.",
            'evaluation': "Correct. This phenomenon is known as the 'leniency bias' and is a well-established finding in educational research. Teachers generally tend to have an optimistic view of their students' abilities, leading to an overestimation of the average student's performance.",
            'is_correct_statement': True
        },
        'B': {
            'statement': "Experienced teacher judgment accuracy correlates with students' achievement around r = .6",
            'evaluation': "Correct. Landmark meta-analyses on this topic (e.g., Hoge & Coladarci, 1989; Südkamp et al., 2012) have found the mean correlation between teacher judgments and standardized achievement tests to be robustly in the range of r = .63 to .66. Therefore, 'around r = .6' is an accurate summary.",
            'is_correct_statement': True
        },
        'C': {
            'statement': "Experienced teachers have more difficulty accurately judging low-performing students.",
            'evaluation': "NOT Correct. This statement is contradicted by research on differential accuracy. Several studies suggest that teachers are often more accurate at identifying students at both the high and low ends of the performance spectrum compared to students in the middle. While there might be a larger overestimation bias for low-performers (level accuracy), their relative rank is often judged accurately (ranking accuracy). Therefore, the blanket statement that they are more difficult to judge is not correct.",
            'is_correct_statement': False
        },
        'D': {
            'statement': "Experienced teachers are typically more accurate in assessing student performance than student teachers.",
            'evaluation': "Correct. Research consistently shows that expertise matters. Experienced teachers tend to have more refined internal standards and more practice, which leads to higher judgment accuracy compared to novice or student teachers.",
            'is_correct_statement': True
        },
        'E': {
            'statement': "Experienced teachers’ judgement accuracy on students achievement is higher than on students motivation.",
            'evaluation': "Correct. Teachers' judgments are generally more accurate for cognitive characteristics like academic achievement, which have direct and observable indicators (e.g., test results, assignments). Judgments of affective or motivational traits rely on interpreting more subtle behaviors and are typically less accurate.",
            'is_correct_statement': True
        }
    }

    incorrect_statement_key = None
    
    print("--- Analysis of Statements on Teacher Judgment Accuracy ---")
    
    for key, value in analysis.items():
        print(f"\nAnalyzing Statement {key}:")
        print(f"Statement: \"{value['statement']}\"")
        # Use textwrap to format the evaluation text nicely
        wrapped_eval = textwrap.fill(f"Evaluation: {value['evaluation']}", width=80)
        print(wrapped_eval)
        if not value['is_correct_statement']:
            incorrect_statement_key = key

    print("\n--- Conclusion ---")
    if incorrect_statement_key:
        print(f"The statement that is NOT correct is {incorrect_statement_key}.")
        print(f"Reasoning: {analysis[incorrect_statement_key]['evaluation']}")
    else:
        print("All statements were found to be correct based on the analysis.")

if __name__ == "__main__":
    analyze_teacher_judgment_statements()
