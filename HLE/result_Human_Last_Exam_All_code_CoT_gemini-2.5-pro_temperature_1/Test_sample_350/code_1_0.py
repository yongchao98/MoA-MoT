import textwrap

def solve_theology_question():
    """
    This function analyzes the provided text about Wolfhart Pannenberg and determines
    the two key historical thinkers he draws upon for his later theology.
    """

    question = """
    Wolfhart Pannenberg... is known, especially in his latter work, for his engage with 
    contemporary physics... he argues for a model of cosmic history that is built upon 
    a contingent concept of time. Upon which two key historic philosophers or 
    theologians does Pannenberg primarily draw upon to support his argument?
    """

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

    # Step 1: Analyze the first key concept: "Theology of History".
    # Pannenberg's entire theological project is a critical dialogue with Georg Hegel's
    # philosophy of history. He adapts Hegel's idea that the meaning of history is
    # revealed at its end. This makes Hegel a definite primary influence.
    # This narrows the possibilities to options E, F, and M.
    think_step1 = "Step 1: Identify the influence for Pannenberg's 'theology of history'."
    conclusion1 = "This is unequivocally Georg Hegel. Pannenberg's work is a direct and critical engagement with Hegel's historical philosophy."
    
    # Step 2: Analyze the second key concept: "Cosmic History" and "Contingency".
    # The prompt mentions Pannenberg's later engagement with physics and a "contingent"
    # concept of time. Hegel's system is often seen as deterministic, not contingent.
    # Friedrich Schelling, however, developed a 'philosophy of nature' and his later work
    # strongly emphasized contingency and freedom, making him the ideal philosophical partner
    # for Pannenberg's project of integrating modern science with theology.
    think_step2 = "Step 2: Identify the influence for 'cosmic history' and 'contingency' to engage with science."
    conclusion2 = "This points directly to Friedrich Schelling and his philosophy of nature."

    # Step 3: Combine the findings to select the correct option.
    think_step3 = "Step 3: Combine the two influences."
    conclusion3 = "The two key figures are Georg Hegel and Friedrich Schelling."
    
    final_choice_letter = 'E'
    final_choice_text = options[final_choice_letter]

    print("Thinking Process:")
    print("-----------------")
    print(textwrap.fill(think_step1 + " " + conclusion1, 80))
    print("\n" + textwrap.fill(think_step2 + " " + conclusion2, 80))
    print("\n" + textwrap.fill(think_step3 + " " + conclusion3, 80))
    print("\n" + "="*25)
    print("Final Answer Selection")
    print("="*25)
    print(f"The correct option is {final_choice_letter}: {final_choice_text}")
    print("The final equation is the combination of these two thinkers:")
    print("Georg Hegel + Friedrich Schelling")


solve_theology_question()
print("<<<E>>>")