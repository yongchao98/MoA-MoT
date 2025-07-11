import itertools

class Parser:
    """
    A recursive descent parser for boolean expressions with precedence rules.
    Grammar:
    E -> T ('|' T)*
    T -> F ('&' F)*
    F -> 'T' | 'F' | '!' F | '(' E ')'
    """
    def __init__(self, text):
        self.text = text
        self.pos = 0

    def parse_expression(self):
        """Parses an expression (handles '|')."""
        node = self.parse_term()
        if node is None:
            return None
        
        while self.pos < len(self.text) and self.text[self.pos] == '|':
            self.pos += 1
            right_node = self.parse_term()
            if right_node is None:
                return None  # Incomplete expression
            node = (node, '|', right_node)
        return node

    def parse_term(self):
        """Parses a term (handles '&')."""
        node = self.parse_factor()
        if node is None:
            return None

        while self.pos < len(self.text) and self.text[self.pos] == '&':
            self.pos += 1
            right_node = self.parse_factor()
            if right_node is None:
                return None # Incomplete expression
            node = (node, '&', right_node)
        return node

    def parse_factor(self):
        """Parses a factor (handles literals, '!', and parentheses)."""
        if self.pos >= len(self.text):
            return None
        
        token = self.text[self.pos]
        
        if token == 'T' or token == 'F':
            self.pos += 1
            return token
        
        if token == '!':
            self.pos += 1
            operand = self.parse_factor()
            if operand is None:
                return None # Incomplete expression
            return ('!', operand)
        
        if token == '(':
            self.pos += 1
            expr = self.parse_expression()
            if expr is None or self.pos >= len(self.text) or self.text[self.pos] != ')':
                return None  # Mismatched or missing parenthesis
            self.pos += 1
            return ('()', expr)
            
        return None

    def get_parse_tree(self):
        """
        Parses the full text. Returns the parse tree if valid and fully consumed,
        otherwise returns None.
        """
        tree = self.parse_expression()
        if tree is not None and self.pos == len(self.text):
            return tree
        return None

def evaluate_tree(tree):
    """Recursively evaluates a parse tree to a boolean value."""
    if tree == 'T':
        return True
    if tree == 'F':
        return False
    
    # Unary operators
    if len(tree) == 2:
        op, operand_tree = tree
        operand_value = evaluate_tree(operand_tree)
        if op == '!':
            return not operand_value
        if op == '()':
            return operand_value
    
    # Binary operators
    if len(tree) == 3:
        left_tree, op, right_tree = tree
        left_value = evaluate_tree(left_tree)
        
        # Short-circuiting for efficiency
        if op == '&' and not left_value:
            return False
        if op == '|' and left_value:
            return True
            
        right_value = evaluate_tree(right_tree)
        
        if op == '&':
            return left_value and right_value
        if op == '|':
            return left_value or right_value
            
    # Should not be reached with a valid tree
    raise ValueError("Invalid parse tree node")

def count_true_expressions():
    """
    Generates all possible strings of length 5, parses them, evaluates them,
    and counts how many are valid and evaluate to True.
    """
    symbols = ['T', 'F', '!', '&', '|', '(', ')']
    true_expression_count = 0
    
    # Generate all 7^5 = 16807 possible strings
    for p in itertools.product(symbols, repeat=5):
        expression_string = "".join(p)
        parser = Parser(expression_string)
        parse_tree = parser.get_parse_tree()
        
        if parse_tree:
            # If the string is a valid expression, evaluate it
            try:
                if evaluate_tree(parse_tree):
                    true_expression_count += 1
            except (ValueError, RecursionError):
                # Catch any evaluation errors from malformed trees, though parser should prevent this
                continue
                
    return true_expression_count

if __name__ == '__main__':
    result = count_true_expressions()
    print(result)

<<<45>>>