import textwrap

def solve_theology_question():
    """
    Analyzes the provided text about Wolfhart Pannenberg and identifies the two key thinkers
    he draws upon for his concept of cosmic history.
    """
    # The text describes Pannenberg's theology of history, his later engagement
    # with physics, and his focus on a "contingent concept of time."

    # 1. Theology of History: Pannenberg's entire theological project is a critical
    #    dialogue with Georg Hegel's philosophy of history. Hegel saw history as the
    #    self-realization of the Absolute Spirit. Pannenberg adopts the universal
    #    scope of history but insists on its contingency and openness to the future,
    #    in contrast to what he saw as Hegel's closed, necessary system.
    thinker1 = "George Hegel"

    # 2. Contingency and Philosophy of Nature: To integrate modern physics and a
    #    "contingent concept of time," Pannenberg turned to another German Idealist,
    #    Friedrich Schelling. Schelling's philosophy of nature (Naturphilosophie)
    #    and his later philosophy of freedom and contingency provided Pannenberg with
    #    the tools to think about God's relationship to a dynamic, evolving, and
    #    non-deterministic cosmos.
    thinker2 = "Friedrich Schelling"

    explanation = (
        "Wolfhart Pannenberg's argument for a 'cosmic history' is built upon a critical "
        "engagement with German Idealism.\n\n"
        "1. For his grand framework of history as the central stage for divine revelation, "
        "he primarily draws upon and critiques {t1}.\n\n"
        "2. For his later work on integrating science and a 'contingent concept of time,' "
        "which emphasizes freedom and nature, he heavily relies on {t2}.\n\n"
        "Therefore, the two key figures are:"
    ).format(t1=thinker1, t2=thinker2)

    print(textwrap.fill(explanation, width=80))
    print("\n---")
    print("Thinker 1: {}".format(thinker1))
    print("Thinker 2: {}".format(thinker2))
    print("---\n")
    print("This corresponds to answer choice E.")

# Execute the function to print the solution.
solve_theology_question()