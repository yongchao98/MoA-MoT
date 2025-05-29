import json

# The solution found by the code
solution = ['fill C', 'pour C->A', 'fill C', 'pour C->A', 'empty A', 'pour C->A', 'fill C', 'pour C->A', 'empty A', 'pour C->A', 'fill C', 'pour C->A', 'empty A', 'pour C->A', 'fill C', 'pour C->A', 'pour C->B', 'pour A->C', 'pour C->B']

# Format the solution as a JSON-parsable list
formatted_solution = json.dumps(solution)
print(f"<<<{formatted_solution}>>>")