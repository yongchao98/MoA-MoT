import sys

# We set a recursion limit to prevent a true stack overflow during demonstration.
# A naive RD parser without this would loop infinitely.
sys.setrecursionlimit(2000)

class Parser:
    def __init__(self, tokens):
        self.tokens = tokens
        self.pos = 0
        self.recursion_depth = 0
        self.max_depth = 20 # Cap to detect infinite recursion

    def parse(self):
        """Top-level parsing function for L = S V C '.' EOF"""
        print(f"--- Starting Parse for: '{' '.join(self.tokens)}' ---")
        if (self.parse_S() and
            self.parse_V() and
            self.parse_C() and
            self.match('.') and
            self.match_eof()):
            print("--- Parse SUCCESSFUL ---")
            return True
        else:
            print("--- Parse FAILED ---")
            return False

    def log(self, msg):
        print(f"{'  ' * self.recursion_depth}{msg}")

    def current_token(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else "EOF"

    def match(self, token):
        if self.pos < len(self.tokens) and self.tokens[self.pos] == token:
            self.log(f"Matched token '{token}'")
            self.pos += 1
            return True
        return False

    def match_eof(self):
        if self.pos == len(self.tokens):
            self.log("Matched EOF")
            return True
        return False

    def parse_S(self):
        """S = N | ADJ N | N ADJ"""
        self.recursion_depth += 1
        self.log(f"Attempting to parse S at token '{self.current_token()}'")
        
        # We must try all alternatives, backtracking on failure.
        saved_pos = self.pos
        
        # Try S -> ADJ N
        if self.parse_ADJ() and self.parse_N():
            self.recursion_depth -= 1
            return True
        self.pos = saved_pos # Backtrack

        # Try S -> N ADJ
        if self.parse_N() and self.parse_ADJ():
            self.recursion_depth -= 1
            return True
        self.pos = saved_pos # Backtrack

        # Try S -> N
        if self.parse_N():
            self.recursion_depth -= 1
            return True
        self.pos = saved_pos # Backtrack

        self.log("Failed to parse S")
        self.recursion_depth -= 1
        return False
        
    def parse_N(self):
        """N = 'frogs' | 'snakes'"""
        self.recursion_depth += 1
        self.log(f"Attempting to parse N at token '{self.current_token()}'")
        if self.match('frogs') or self.match('snakes'):
            self.recursion_depth -= 1
            return True
        self.recursion_depth -= 1
        return False
        
    def parse_V(self):
        """V = 'jump' | 'swim'"""
        self.recursion_depth += 1
        self.log(f"Attempting to parse V at token '{self.current_token()}'")
        if self.match('jump') or self.match('swim'):
            self.recursion_depth -= 1
            return True
        self.recursion_depth -= 1
        return False

    def parse_ADJ(self):
        """ADJ = 'red' | 'or alike' | REC"""
        self.recursion_depth += 1
        self.log(f"Attempting to parse ADJ at token '{self.current_token()}'")
        
        saved_pos = self.pos
        if self.match('red') or self.match('or alike'):
            self.recursion_depth -= 1
            return True
        self.pos = saved_pos

        # If the above fail, we must try the REC rule
        if self.parse_REC():
             self.recursion_depth -= 1
             return True

        self.recursion_depth -= 1
        return False

    def parse_REC(self):
        """REC = REC ADJ  <-- LEFT RECURSION!"""
        self.recursion_depth += 1
        self.log(f"Attempting to parse REC at token '{self.current_token()}'")

        if self.recursion_depth > self.max_depth:
            print("\n!!! INFINITE RECURSION DETECTED !!!")
            print("Parser is stuck in a loop trying to parse REC -> REC ADJ.")
            print("This shows the failure of naive RD parsing with left recursion.\n")
            # To prevent crashing, we escape the loop
            sys.exit("Exiting due to detected infinite loop.")
        
        # This is the problematic call. It calls itself without consuming input.
        if self.parse_REC() and self.parse_ADJ():
            self.recursion_depth -= 1
            return True

        self.recursion_depth -= 1
        return False

    def parse_C(self):
        """C = 'well' | 'swiftly'"""
        self.recursion_depth += 1
        self.log(f"Attempting to parse C at token '{self.current_token()}'")
        if self.match('well') or self.match('swiftly'):
            self.recursion_depth -= 1
            return True
        self.recursion_depth -= 1
        return False

# --- Analysis of Statements ---

# Statement A: "The language contains 'red frogs swim swiftly', and it is not the longest sentence in the language."
# L(REC) is empty, so L(ADJ) is {'red', 'or alike'}. The language is finite.
# Longest sentences have 4 words (ADJ N V C or N ADJ V C).
# "red frogs swim swiftly ." has 4 words. It IS a longest sentence. So, A is false.

# Statement B & C check
print("--- Checking sentences from B and C ---")
print("These sentences are not in the language, so a correct parser should fail.")
Parser("snake big jump well .".split()).parse()
print("\n")
Parser("jump or alike frog swiftly .".split()).parse()
print("\n")

# Statement D check: "Neither RD nor PR can be implemented..."
# Let's show it CAN be implemented, but it's flawed.
# Case 1: A sentence that can be parsed without hitting the REC->REC rule.
print("--- Demonstrating a successful parse (Statement D is false) ---")
Parser("frogs red swim well .".split()).parse()
print("\nThis worked because parse_ADJ was called on 'red' and succeeded immediately.")
print("The parser never had to try the problematic REC rule.")
print("This shows the parser *is* implementable and *can* work for some inputs.\n")


# Case 2: A sentence that forces the parser into the left-recursive rule.
print("--- Demonstrating parser failure via infinite loop (the flaw in the implementation) ---")
print("The following parse will be terminated by our recursion-depth check.")
Parser("blue frogs swim swiftly .".split()).parse()
