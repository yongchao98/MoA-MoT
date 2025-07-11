import textwrap

def solve_pannenberg_question():
    """
    This function analyzes the theological question about Wolfhart Pannenberg
    and identifies the primary influences on his later thought.
    """

    # A dictionary mapping the answer choices to their content
    options = {
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

    # The correct answer key based on theological scholarship.
    correct_key = 'E'

    # The explanation for the answer.
    explanation = """
    The question asks about the primary influences on Wolfhart Pannenberg's concept of a 'cosmic history' that engages with modern physics. The analysis points to two key figures from German Idealism:

    1.  Georg Hegel: Pannenberg's framework of understanding history as the primary mode of God's universal self-revelation is a critical inheritance from Hegel's philosophy of history.

    2.  Friedrich Schelling: To integrate natural science and physics into this historical framework, Pannenberg drew heavily on Schelling's 'Naturphilosophie' (philosophy of nature). Schelling viewed nature itself as a dynamic, historical process, which provided Pannenberg with the philosophical tools to speak of a single 'cosmic history' encompassing both human and natural events.
    
    Therefore, the combination of Hegel and Schelling is the most accurate answer.
    """

    print("Analysis and Solution:")
    # Use textwrap to format the explanation nicely.
    print(textwrap.dedent(explanation).strip())
    print("\n" + "="*40)
    print("Conclusion:")
    print(f"The correct option is {correct_key}.")
    print(f"Selected Answer: {options[correct_key]}")
    print("="*40)


solve_pannenberg_question()