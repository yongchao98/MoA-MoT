def main_solution(n):
    import copy
    
    var_list = []
    statements = []
    
    for i in range(n):
        var_list.append(["a{}".format(i), "~a{}".format(i)])
    
    helper(var_list, 0, statements, [])
    
    fun = lambda x: " & ".join(x)
    fun2 = lambda x: "({})".format(x)
    statements = map(fun, statements)
    statements = map(fun2, statements)
    return " | ".join(statements)


def helper(v_list, i, statements, statement):
    if i == len(v_list):
        t = copy.copy(statement)
        statements.append(t)
        return
    
    statement.append(v_list[i][0])
    helper(v_list, i+1, statements, statement)
    statement.pop()
    statement.append(v_list[i][1])
    helper(v_list, i+1, statements, statement)
    statement.pop()

# Generate the DNF for n = 5
print(main_solution(5))