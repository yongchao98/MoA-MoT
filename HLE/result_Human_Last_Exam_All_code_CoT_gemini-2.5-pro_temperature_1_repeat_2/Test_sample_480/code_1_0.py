import sys

def solve_natural_logic_inference():
    """
    Solves the entailment problem using MacCartney's natural logic framework.
    """
    # Define the 7 semantic relations and their names
    relation_names = {
        '=': 'equivalence',
        'sq': 'forward entailment (subsumption)',
        'sp': 'reverse entailment (subsumption)',
        '^': 'negation (contradiction)',
        '|': 'alternation (mutual exclusion)',
        '~': 'cover (exhaustion)',
        '#': 'independence',
    }

    # MacCartney's join table for composing relations (R1 o R2)
    # This is join(R_old, R_new_projected)
    join_table = {
        # R_new ->
        # R_old v   =    sq   sp   ^    |    ~    #
        '=':      {'=': '=', 'sq': 'sq', 'sp': 'sp', '^': '^', '|': '|', '~': '~', '#': '#'},
        'sq':     {'=': 'sq', 'sq': 'sq', 'sp': '#', '^': '|', '|': '|', '~': '#', '#': '#'},
        'sp':     {'=': 'sp', 'sq': '#', 'sp': 'sp', '^': '~', '|': '#', '~': '~', '#': '#'},
        '^':      {'=': '^', 'sq': '~', 'sp': '|', '^': '=', '|': 'sp', '~': 'sq', '#': '#'},
        '|':      {'=': '|', 'sq': '#', 'sp': '|', '^': 'sq', '|': '#', '~': 'sq', '#': '#'},
        '~':      {'=': '~', 'sq': '~', 'sp': '#', '^': 'sp', '|': 'sp', '~': '#', '#': '#'},
        '#':      {'=': '#', 'sq': '#', 'sp': '#', '^': '#', '|': '#', '~': '#', '#': '#'}
    }

    # Projection table for lexical relations in different contexts.
    # We only need UP ('+') and DOWN ('-') contexts for this problem.
    # DOWN corresponds to the 'not' column in MacCartney's table.
    projection_table = {
        # Lexical v  Context ->
        '=':      {'+': '=', '-': '='},
        'sq':     {'+': 'sq', '-': 'sp'},
        'sp':     {'+': 'sp', '-': 'sq'},
        '^':      {'+': '^', '-': '^'},
        '|':      {'+': '|', '-': '~'},
        '~':      {'+': '~', '-': '|'},
        '#':      {'+': '#', '-': '#'},
    }

    # Sequence of edits to transform Premise to Hypothesis, ordered left-to-right
    edits = [
        {'desc': "'is singing' -> 'is not singing'", 'lexical': '^'},
        {'desc': "'a pop song' -> 'a song'", 'lexical': 'sq'},
        {'desc': "'Taylor Swift' -> 'Michael Jackson'", 'lexical': '|'},
    ]

    print("Premise: 'Mark is singing a pop song by Taylor Swift'")
    print("Hypothesis: 'Mark is not singing a song by Michael Jackson'\n")
    print("Calculating final entailment relation step-by-step:\n")

    # Initial state
    current_relation = '='
    context = '+' # Initial context is upward monotone

    # Process each edit
    for i, edit in enumerate(edits):
        print(f"Step {i+1}: Edit {edit['desc']}")
        
        lexical_relation = edit['lexical']
        print(f"  - Lexical relation: '{lexical_relation}' ({relation_names[lexical_relation]})")

        context_str = 'UP' if context == '+' else 'DOWN'
        print(f"  - Current context: {context_str}")

        projected_relation = projection_table[lexical_relation][context]
        print(f"  - Projected relation: '{projected_relation}' ({relation_names[projected_relation]})")

        old_relation = current_relation
        current_relation = join_table[old_relation][projected_relation]
        
        print(f"  - Composed relation: join('{old_relation}', '{projected_relation}') = '{current_relation}' ({relation_names[current_relation]})")
        print("-" * 20)

        # The introduction of negation makes all subsequent contexts downward monotone
        if lexical_relation == '^':
            context = '-'

    final_relation_name = relation_names[current_relation]
    print(f"\nThe final projected natural logic operator for the inference is '{current_relation}'.")
    print(f"The name of this operator is: {final_relation_name.split(' (')[0]}")
    
    # Set the return value for the special answer format
    # This part is for the platform, the printed output is for the user.
    sys.stdout.final_answer = final_relation_name.split(' (')[0]


# Run the simulation
# In a real scenario, this would be the main execution block.
# We will capture its output.
class Capturing:
    def __init__(self):
        self._stdout = None
        self.final_answer = None

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._stdout

    def write(self, data):
        self._stdout.write(data)

with Capturing() as output:
    solve_natural_logic_inference()

# The final answer is extracted from the function's side effect.
final_answer = output.final_answer
print(f'<<<{final_answer}>>>')
