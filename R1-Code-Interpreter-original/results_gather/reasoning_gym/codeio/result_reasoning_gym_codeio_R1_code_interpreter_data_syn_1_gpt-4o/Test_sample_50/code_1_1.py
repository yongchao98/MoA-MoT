def subset_sum(conjunto, capacidade):
    n = len(conjunto)
    dp = [[False] * (capacidade + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = True

    for i in range(1, n + 1):
        for j in range(1, capacidade + 1):
            if conjunto[i - 1] > j:
                dp[i][j] = dp[i - 1][j]
            else:
                dp[i][j] = dp[i - 1][j] or dp[i - 1][j - conjunto[i - 1]]

    return dp

def find_subsets(conjunto, dp, capacidade):
    def backtrack(i, j, path):
        if j == 0:
            result.append(path)
            return
        if i == 0:
            return
        if dp[i - 1][j]:
            backtrack(i - 1, j, path)
        if j >= conjunto[i - 1] and dp[i - 1][j - conjunto[i - 1]]:
            backtrack(i - 1, j - conjunto[i - 1], path + [conjunto[i - 1]])

    result = []
    backtrack(len(conjunto), capacidade, [])
    return result

def main_solution(conjunto, capacidade, tecninca, todos):
    if tecninca == 1:
        dp = subset_sum(conjunto, capacidade)
        if dp[len(conjunto)][capacidade]:
            subsets = find_subsets(conjunto, dp, capacidade)
            return subsets if todos else subsets[:1]
    return []

# Given input
input_data = {'conjunto': [17, 10, 6, 6, 9, 6, 2, 14], 'capacidade': 61, 'tecninca': 1, 'todos': True}
result = main_solution(input_data['conjunto'], input_data['capacidade'], input_data['tecninca'], input_data['todos'])
print(result)