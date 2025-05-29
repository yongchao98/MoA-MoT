# Define sets for butterflies, animals, and lions
butterflies = set(['b1', 'b2', 'b3'])  # Example butterflies
animals = set(['b1', 'b2'])  # Some butterflies are animals
lions = set()  # No animals are lions

# Check if some butterflies are not lions
some_butterflies_not_lions = any(butterfly not in lions for butterfly in butterflies)

print(some_butterflies_not_lions)