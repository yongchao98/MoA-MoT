# Define the necessary classes and functions
class Stream(object):
    def __init__(self, string, position=0):
        self.string = string
        self.position = position
        
    def __repr__(self):
        return "'%s'@%i" % (self.string, self.position)
    
    def advanced(self, characters):
        return Stream(self.string, self.position + characters)
    
    def startswith(self, prefix):
        return self.string[self.position:].startswith(prefix)

class Grammer(object):
    def __getitem__(self, name):
        return self.__getattribute__(name)
    
    def start(self):
        raise Exception("Start must be implemented to create a grammer.")

    def parse(self, s):
        if not isinstance(s, Stream):
            s = Stream(s)
        res = self.start().apply(s, self)
        return res

class Parser(object):
    def __init__(self):
        self.postprocessor = lambda x: x
    
    def apply(self, stream, grammer):
        raise Exception("apply called on generic Parser.")

class Word(Parser):
    def __init__(self, symbol, post=None):
        super(Word, self).__init__()
        self.symbol = symbol
        if post:
            self.postprocessor = post
    
    def __repr__(self):
        return "'%s'" % self.symbol
    
    def apply(self, stream, grammer):
        if stream.startswith(self.symbol):
            return (stream.advanced(len(self.symbol)), True, self.postprocessor(None))
        else:
            return (stream, False, None)

class Charset(Parser):
    def __init__(self, charset):
        super(Charset, self).__init__()
        self.charset = charset
    
    def __repr__(self):
        return "[%s]" % str(self.charset)
    
    def apply(self, stream, grammer):
        adv = 0
        rest = stream.string[stream.position:]
        while adv < len(rest) and rest[adv] in self.charset:
            adv += 1
        if adv > 0:
            ret = rest[:adv]
            return (stream.advanced(adv), True, self.postprocessor(ret))
        else:
            return (stream, False, None)

class AlphaSet(object):
    def __repr__(self):
        return "a-zA-z"
    
    def __contains__(self, char):
        return ("a" <= char <= "z") or ("A" <= char <= "Z")
    
def Alpha():
    return Charset(AlphaSet())

# Main function to solve the problem
def main_solution(input_string, parser_name):
    # Define a simple grammar for parsing
    class SimpleGrammar(Grammer):
        def __init__(self):
            self.word = Word("hello")
            self.alpha = Alpha()
    
    # Create an instance of the grammar
    grammar = SimpleGrammar()
    
    # Convert input string to Stream
    stream = Stream(input_string)
    
    # Apply the parser specified by parser_name
    parser = getattr(grammar, parser_name)
    new_stream, matched, production = parser.apply(stream, grammar)
    
    # Return the result as a dictionary
    return {
        "new_stream": repr(new_stream),
        "matched": matched,
        "production": production
    }

# Execute the main function with the given input
result = main_solution('hellozO29825', 'alpha')
print(result)