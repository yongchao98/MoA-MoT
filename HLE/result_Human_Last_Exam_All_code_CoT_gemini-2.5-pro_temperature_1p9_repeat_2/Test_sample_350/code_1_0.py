import sys

def solve_pannenberg_question():
    """
    Analyzes the influences on Wolfhart Pannenberg's theology of history
    and identifies the correct philosophers from the given choices.
    """
    question = """
    Wolfhart Pannenberg... argues for a model of cosmic history that is built
    upon a contingent concept of time. Upon which two key historic philosophers
    or theologians does Pannenberg primarily draw upon to support his argument?
    """

    # Analysis: Pannenberg's theology of history is famously a critical re-engagement
    # with German Idealism. His idea of history as God's self-revelation is a
    # direct response to Hegel. His later work integrating science, physics, and
    # contingency draws heavily on Schelling's philosophy of nature (Naturphilosophie).
    # Therefore, Hegel and Schelling are the primary influences for this specific argument.

    answer_choices = {
        'A': 'Dietrich Bonhoeffer and Jurgen Moltmann',
        'B': 'Paul Tillich and Karl Barth',
        'C': 'Martin Luther and Friedrich Schelling',
        'D': 'Martin Luther and Rudolf Bultmann',
        'E': 'George Hegel and Friedrich Schelling',
        'F': 'Georg Hegel and Martin Heidegger',
        'G': 'John Duns Scotus and Rudolf Bultmann',
        'H': 'Gottfried Leibniz and Martin Heidegger',
        'I': 'Paul Tillich and Thomas Aquinas',
        'J': 'Martin Luther and Martin Heidegger',
        'K': 'Paul Tillich and Jurgen Moltmann',
        'L': 'Gottfried Leibniz and Jurgen Moltmann',
        'M': 'George Hegel and Gottfried Leibniz',
        'N': 'John Duns Scotus and Paul Tillich',
        'O': 'Paul Tillich and Friedrich Schelling',
        'P': 'Dietrich Bonhoeffer and Rudolf Bultmann',
        'Q': 'Martin Luther and Jurgen Moltmann',
        'R': 'Karl Barth and Friedrich Schelling',
        'S': 'Paul Tillich and John Calvin',
        'T': 'John Duns Scotus and Friedrich Schelling',
        'U': 'Paul Tillich and Martin Heidegger',
        'V': 'Martin Luther and John Calvin',
        'W': 'Dietrich Bonhoeffer and Thomas Aquinas',
        'X': 'Karl Barth and Rudolf Bultmann',
        'Y': 'John Duns Scotus and Martin Heidegger',
        'Z': 'Martin Luther and Thomas Aquinas'
    }

    correct_answer_key = 'E'
    influencer_1 = "Georg Hegel"
    influencer_2 = "Friedrich Schelling"

    print(f"The two key thinkers Pannenberg draws upon are: {influencer_1} and {influencer_2}.")
    print(f"This corresponds to answer choice: {correct_answer_key}")
    print(f"Final Answer: {correct_answer_key}. {answer_choices[correct_answer_key]}")


solve_pannenberg_question()