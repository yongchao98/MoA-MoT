def print_solution(moves):
    result = "<<<" + "\n".join(moves) + ">>>"
    print(result)

moves = [
    "Move D from 1 to 4",
    "Move I from 1 to 3",
    "Move G from 1 to 4",
    "Move H from 1 to 3",
    "Move B from 1 to 2",
    "Move F from 2 to 3",
    "Move D from 4 to 2",
    "Move G from 4 to 2",
    "Move C from 4 to 3",
    "Move J from 4 to 3",
    "Move J from 3 to 4",
    "Move C from 3 to 4",
    "Move F from 3 to 4",
    "Move I from 3 to 4",
    "Move H from 3 to 1",
    "Move I from 4 to 3",
    "Move F from 4 to 3",
    "Move C from 4 to 3",
    "Move J from 4 to 3",
    "Move G from 2 to 4",
    "Move D from 2 to 4",
    "Move G from 4 to 2",
    "Move D from 4 to 2"
]

print_solution(moves)