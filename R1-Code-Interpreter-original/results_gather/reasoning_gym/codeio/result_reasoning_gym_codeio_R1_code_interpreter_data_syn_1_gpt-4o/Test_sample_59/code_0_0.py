def is_pareto_improvement(agents_values, option1, option2):
    improvement = False
    for agent in agents_values:
        if agent[option1] < agent[option2]:
            return False
        if agent[option1] > agent[option2]:
            improvement = True
    return improvement

def is_pareto_optimal(agents_values, option):
    num_options = len(agents_values[0])
    for other_option in range(num_options):
        if other_option == option:
            continue
        if all(agent[option] >= agent[other_option] for agent in agents_values) and any(agent[option] > agent[other_option] for agent in agents_values):
            return False
    return True

def main_solution(agents_values, option1, option2, option):
    is_improvement = is_pareto_improvement(agents_values, option1, option2)
    is_optimal = is_pareto_optimal(agents_values, option)
    return {"is_improvement": is_improvement, "is_optimal": is_optimal}

# Given input
agents_values = [
    [40.22501628786998, 45.898115518933544, 78.2210370468844, 44.01180108297291],
    [45.807377865170174, 25.555509510392625, 18.789633067917432, 84.03256563139756],
    [80.0154235142718, 21.77527842219694, 79.1500808422558, 38.64861186948]
]
option1 = 0
option2 = 2
option = 3

# Calculate the result
result = main_solution(agents_values, option1, option2, option)
print(result)