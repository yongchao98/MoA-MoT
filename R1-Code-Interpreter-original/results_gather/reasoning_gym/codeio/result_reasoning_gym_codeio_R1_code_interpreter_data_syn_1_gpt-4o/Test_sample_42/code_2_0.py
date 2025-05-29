def main_solution(n):
    from itertools import product
    
    var_list = ["a{}".format(i) for i in range(n)]
    statements = []
    
    for combination in product([True, False], repeat=n):
        term = []
        for i, value in enumerate(combination):
            if value:
                term.append(var_list[i])
            else:
                term.append("~" + var_list[i])
        statements.append(" & ".join(term))
    
    return " | ".join(f"({s})" for s in statements)

# Generate the DNF for n = 5
print(main_solution(5))