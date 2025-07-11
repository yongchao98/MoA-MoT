import textwrap

def explain_bacterial_resistance():
    """
    Explains the reasoning behind the correct answer to the biology question.
    """
    explanation = """
    The central problem is to explain how a bacterium with a stable genome (Bacterium 2), which evolves mainly through mutation, can acquire drug resistance as quickly as a bacterium that uses rapid lateral gene transfer (LGT, Bacterium 1).

    1.  **The Baseline (Bacterium 1):** LGT is very fast. It allows for the transfer of pre-existing resistance genes (often on plasmids) to many cells in a population simultaneously. This is the benchmark for "fast" resistance acquisition.

    2.  **The Challenge (Bacterium 2):** Evolution by mutation is a multi-step process: a rare mutation must occur, and then the resulting mutant must spread through the population via natural selection. For this to be as fast as LGT, the spread of the mutant must be extremely rapid.

    3.  **Evaluating the Options:**
        *   **A is insufficient:** A rare mutation alone doesn't explain the rapid *spread*.
        *   **C is a distraction:** It suggests experimental error, not a biological mechanism.
        *   **D is unlikely:** Without compensatory mutations, the initial resistance mutation often carries a fitness cost, which would *slow down* its spread, not speed it up.
        *   **E is correct but incomplete:** It mentions compensatory mutations but misses other key accelerators.

    4.  **Why B is the Best Answer:** Option B provides the most complete and powerful explanation. Let's break it down:
        *   **'rare resistance mutations':** This is the necessary starting point.
        *   **'compensatory mutations that increased the fitness to a great extent':** This is the crucial part. Resistance mutations can be costly. Compensatory mutations eliminate that cost, making the resistant bacterium highly fit and able to rapidly outcompete and replace the original, susceptible population. This explains the *speed*.
        *   **'also led to cross-resistance':** This is a significant accelerator. If the same mutations confer resistance to multiple drugs, the bacterium gains a massive selective advantage in a multi-drug environment. This high-fitness, multi-resistant strain can sweep through the population at a pace that rivals the spread of a resistance plasmid via LGT.

    Therefore, the combination of highly effective compensatory mutations and the bonus of cross-resistance provides a robust mechanism for the observed rapid evolution in the second bacterium.
    """
    print(textwrap.dedent(explanation).strip())

explain_bacterial_resistance()