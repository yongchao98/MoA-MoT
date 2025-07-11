import textwrap

def analyze_ontario_non_compete_law():
    """
    Analyzes which employee can have a valid non-compete clause in Ontario as of Jan 2023.
    """

    # Step 1: Explain the relevant law in Ontario.
    print("Step 1: Understanding Ontario's Law on Non-Competition Clauses")
    explanation = """
    In Ontario, the Employment Standards Act, 2000 (ESA) was amended by the Working for Workers Act, 2021.
    Effective October 25, 2021, employers are prohibited from entering into employment agreements that include non-competition clauses.
    Any non-compete clause signed on or after this date is generally void.
    """
    print(textwrap.dedent(explanation))

    # Step 2: Explain the key exception to the rule.
    print("Step 2: Identifying the Exception to the Rule")
    exception_explanation = """
    There is a significant exception to this ban. The prohibition does NOT apply to employees who are "executives".
    The ESA narrowly defines an "executive" as someone holding a C-suite title (e.g., Chief Executive Officer, President, Chief Financial Officer) or "any other chief executive position."
    This means the exception is intended for the highest level of company leadership.
    """
    print(textwrap.dedent(exception_explanation))

    # Step 3: Analyze each employee option against the "executive" exception.
    print("Step 3: Analyzing the Employee Options")
    employees = {
        "A": "A Restaurant Manager at a local restaurant. This is a management role, but it does not meet the high threshold of a C-suite level 'executive' for the entire business. The non-compete clause would be void.",
        "B": "An Associate Lawyer at a large corporate law firm. This is a professional role, not an executive one. Furthermore, the Law Society of Ontario's Rules of Professional Conduct also prohibit non-competes for lawyers. The clause would be void.",
        "C": "A Branch Manager at a Schedule 1 bank. This is the most senior position among the choices. While a Branch Manager is not a C-suite executive of the entire bank, they hold the top leadership position for their specific branch. In the context of the choices provided, this is the only role that could plausibly be argued to be a 'chief executive position' of a business unit, thus potentially fitting the exception.",
        "D": "A Hairdresser at a local salon. This is a skilled service provider, not an executive. The non-compete clause would be void.",
        "E": "A Cashier at a large franchise restaurant. This is a front-line, non-management role and is clearly not an executive. The non-compete clause would be void."
    }

    for key, value in employees.items():
        print(f"\nAnalysis of Option {key}:")
        print(textwrap.fill(value, width=80))

    # Step 4: Conclude with the most likely answer.
    print("\n" + "="*40)
    print("Conclusion")
    print("="*40)
    conclusion = """
    Based on the analysis, the roles of Restaurant Manager, Associate Lawyer, Hairdresser, and Cashier do not meet the narrow definition of "executive" under the ESA. The Branch Manager at a major bank is the most senior role and the only one that could potentially be argued to fall under the "executive" exemption, making it the only viable answer among the choices.
    """
    print(textwrap.dedent(conclusion))
    
    final_answer = "C"
    print(f"The final answer is {final_answer}")


if __name__ == "__main__":
    analyze_ontario_non_compete_law()
    print("<<<C>>>")