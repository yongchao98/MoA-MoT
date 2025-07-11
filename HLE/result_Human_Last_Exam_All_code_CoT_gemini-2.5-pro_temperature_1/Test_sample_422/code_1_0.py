import sys

# Set a lower recursion limit to prevent a long crash when demonstrating the regress problem.
# The default is often 1000.
sys.setrecursionlimit(15)

class Belief:
    """A simple class to represent a belief and its properties."""
    def __init__(self, content, is_true, justification=None):
        self.content = content
        self.is_true = is_true
        # The justification for a belief is often composed of other beliefs.
        self.justification = justification

def run_analysis():
    """Demonstrates and explains the two problems with the JTB definition."""
    print("Analyzing the JTB definition: Knowledge = Justified True Belief")
    print("Constraint: The only available epistemic states are 'Knowledge' and 'Belief'.")
    print("="*60)

    # --- Problem 1: The Infinite Regress of Justification ---
    print("Problem 1: The Infinite Regress of Justification\n")
    print("The JTB model requires justification. But what is the status of the justification itself?")
    print("If we argue that the justification must also be 'Known', we run into an infinite regress.")
    print("Consider this function:")
    print("""
    def is_knowledge_regress(belief):
        # A belief is known only if it's true and its justification is also known.
        if not belief.is_true:
            return False
        if belief.justification is None: # Unjustified
            return False
        # Recursive step: Check if the justification is known
        return is_knowledge_regress(belief.justification)
    """)

    # Let's create a chain of beliefs to demonstrate this.
    j2 = Belief(content="Evidence E exists", is_true=True, justification=None) # The chain must end somewhere
    j1 = Belief(content="Proposition Q is true", is_true=True, justification=j2)
    p = Belief(content="Proposition P is true", is_true=True, justification=j1)
    
    # We add a base case to our real demonstration function to avoid a crash,
    # but highlight that this base case itself is problematic.
    def demonstrate_regress(belief, depth=0):
        if belief.justification is None:
            print(f"{'  '*depth}--> Reached a foundational belief ('{belief.content}'). But is this foundation known or just assumed? This is the problem.")
            return False # An unjustified belief cannot be knowledge.
        print(f"{'  '*depth}Checking belief: '{belief.content}'. To know this, we must know its justification...")
        return demonstrate_regress(belief.justification, depth + 1)
        
    print("Executing this logic shows that to know P, you must know J1; to know J1, you must know J2, and so on.")
    demonstrate_regress(p)
    print("\nThe chain can either go on forever (infinite regress) or stop at an unjustified belief, meaning no knowledge could ever be truly established.\n")
    
    print("="*60)

    # --- Problem 2: The Devaluation of Justified False Beliefs ---
    print("Problem 2: The Devaluation of Justified False Beliefs\n")
    print("The JTB model fails to distinguish between a well-reasoned but false belief and a mere guess.")
    print("If a belief is not Knowledge (not a JTB), it falls into the only other category: 'Belief'.")
    
    # A simple checker for this demonstration
    def is_knowledge_simple(belief):
        return belief.is_true and (belief.justification is not None)

    # Case A: A belief that is justified, but happens to be false.
    # e.g., Believing the Ptolemaic model based on the extensive evidence and complex math available for 1500 years.
    justified_false_belief = Belief(
        content="The sun revolves around the Earth.",
        is_true=False,
        justification="Based on the highly-developed Ptolemaic system and naked-eye observations."
    )

    # Case B: An unjustified belief (a wild guess).
    unjustified_belief = Belief(
        content="The moon is made of green cheese.",
        is_true=False,
        justification=None
    )

    is_knowledge_A = is_knowledge_simple(justified_false_belief)
    is_knowledge_B = is_knowledge_simple(unjustified_belief)

    print(f"Case A: A Justified but False Belief")
    print(f"   - Content: '{justified_false_belief.content}'")
    print(f"   - Is it Knowledge? -> {is_knowledge_A}")

    print(f"\nCase B: An Unjustified Belief (a guess)")
    print(f"   - Content: '{unjustified_belief.content}'")
    print(f"   - Is it Knowledge? -> {is_knowledge_B}")

    print("\nResult: Both are classified as 'not Knowledge'.")
    print("By forcing both into the single category of 'Belief', the model fails to recognize the different epistemic value between a responsibly held belief (A) and a wild guess (B).")

if __name__ == '__main__':
    run_analysis()