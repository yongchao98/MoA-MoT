import inspect

class Object:
    """Represents a singular object/term, like 'a'."""
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return self.name

class Predicate:
    """Represents a predicate concept, like 'F'."""
    def __init__(self, symbol):
        self.symbol = symbol
    def __str__(self):
        return self.symbol

class UniversalQuantifier:
    """Represents the concept of universal quantification, '∀'."""
    def __str__(self):
        return "∀"

class Proposition:
    """Represents a simple proposition, e.g., F(a)."""
    def __init__(self, predicate, obj):
        if not isinstance(predicate, Predicate) or not isinstance(obj, Object):
            raise TypeError("Proposition must be formed by a Predicate and an Object.")
        self.predicate = predicate
        self.obj = obj
        # As per the Generality Constraint, understanding this proposition
        # means we have access to its constituent parts.
        self.conceptual_parts = {'predicate': self.predicate, 'object': self.obj}

    def __str__(self):
        return f"{self.predicate}({self.obj})"

class QuantifiedProposition:
    """Represents a universally quantified proposition, e.g., ∀x(F(x))."""
    def __init__(self, quantifier, predicate, variable='x'):
        if not isinstance(quantifier, UniversalQuantifier) or not isinstance(predicate, Predicate):
            raise TypeError("QuantifiedProposition needs a Quantifier and a Predicate.")
        self.quantifier = quantifier
        self.predicate = predicate
        self.variable = variable
        self.conceptual_parts = {'quantifier': self.quantifier, 'predicate': self.predicate, 'variable': self.variable}


    def __str__(self):
        return f"{self.quantifier}{self.variable}({self.predicate}({self.variable}))"

# --- Main Demonstration ---

# 1. Start with the premises.
# Premise A: We understand the proposition "Fa". Let's create it.
# This implies we possess the concepts 'F' and 'a'.
predicate_F = Predicate("F")
object_a = Object("a")
prop_Fa = Proposition(predicate_F, object_a)
print(f"Assume we understand the initial proposition: {prop_Fa}")

# Premise B: We independently understand universal quantification.
quantifier_forall = UniversalQuantifier()
print(f"Assume we also understand the concept of the universal quantifier: {quantifier_forall}")
print("-" * 20)

# 2. Apply the Generality Constraint.
# Understanding 'Fa' means we can abstract its parts. Let's get the predicate.
abstracted_predicate = prop_Fa.conceptual_parts['predicate']
print(f"From '{prop_Fa}', we can access the conceptual component part: Predicate '{abstracted_predicate}'")

# 3. Recombine concepts to form a new thought.
# We now combine the abstracted predicate with our concept of the quantifier.
new_quantified_prop = QuantifiedProposition(quantifier_forall, abstracted_predicate)
print(f"By recombining '{abstracted_predicate}' with '{quantifier_forall}', we form the new proposition: {new_quantified_prop}")
print("-" * 20)

# 4. As requested, output each component "number" of the final "equation".
final_prop = new_quantified_prop
print(f"The final proposition is: {final_prop}")
print("Its conceptual components are:")
for part_name, part_value in final_prop.conceptual_parts.items():
    print(f"- {part_name.capitalize()}: {part_value}")
