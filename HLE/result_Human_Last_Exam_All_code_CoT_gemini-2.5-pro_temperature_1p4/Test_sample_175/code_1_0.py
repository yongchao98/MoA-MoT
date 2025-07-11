import textwrap

def analyze_automation_bias():
    """
    Analyzes options related to mitigating automation bias in education and identifies
    the one that would NOT be an effective measure.
    """
    print("Analyzing which measure will NOT reduce automation bias in teachers' assessments of student performance:")
    print("-" * 80)

    # Option A Analysis
    option_a_text = """A. Encouraging teachers accountability for decisions made with AI support.
    Analysis: If teachers are held accountable for the final assessment, they are motivated to critically evaluate the AI's suggestion rather than accepting it passively. This increases their engagement and vigilance.
    Effect: Reduces automation bias."""
    print(textwrap.dedent(option_a_text))
    print("-" * 80)

    # Option B Analysis
    option_b_text = """B. Displays the confidence level of the AI update after each decision versus the overall confidence level of the AI.
    Analysis: Providing the AI's confidence level gives teachers more information to judge the reliability of the output. A low confidence score would signal that the AI's suggestion should be carefully scrutinized.
    Effect: Reduces automation bias."""
    print(textwrap.dedent(option_b_text))
    print("-" * 80)

    # Option C Analysis
    option_c_text = """C. Regular practice using AI tools to assess student performance.
    Analysis: While practice can sometimes lead to complacency, properly designed training would expose teachers to the AI's strengths and weaknesses, including its errors. This helps them learn when to trust the AI and when to be skeptical. It is generally considered a way to improve user performance with a system.
    Effect: Aims to reduce automation bias (by improving calibration of trust)."""
    print(textwrap.dedent(option_c_text))
    print("-" * 80)

    # Option D Analysis
    option_d_text = """D. Making the AI advice more salient on the interface.
    Analysis: Salience means making something more prominent or noticeable (e.g., using a larger font, bright colors, or a central position). Making the AI's advice more salient gives it greater psychological weight, making the teacher more likely to notice and accept it. This encourages reliance on the AI as a shortcut, which is the definition of automation bias.
    Effect: Increases (or does not reduce) automation bias."""
    print(textwrap.dedent(option_d_text))
    print("-" * 80)

    # Option E Analysis
    option_e_text = """E. Requiring teachers to justify decisions made based on AI suggestions.
    Analysis: Similar to accountability, requiring justification forces teachers to engage in deeper cognitive processing. They cannot simply click 'agree' but must formulate a reason for their decision, which involves actively thinking about the problem.
    Effect: Reduces automation bias."""
    print(textwrap.dedent(option_e_text))
    print("-" * 80)

    print("\nConclusion: Making the AI advice more salient (Option D) would encourage over-reliance, not reduce it. All other options are strategies designed to promote critical thinking and user vigilance.")
    print("<<<D>>>")

analyze_automation_bias()