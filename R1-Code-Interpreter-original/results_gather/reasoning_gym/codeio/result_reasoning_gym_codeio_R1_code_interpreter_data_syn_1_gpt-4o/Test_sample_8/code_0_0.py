from collections import defaultdict

def get_parent(a, parent):
    if a != parent[a]:
        parent[a] = get_parent(parent[a], parent)
    return parent[a]

def set_parent(a, b, parent):
    p1 = get_parent(a, parent)
    p2 = get_parent(b, parent)
    if p1 != p2:
        parent[p1] = p2

def is_bipartite(parent, adj, n):
    color = {}
    for i in range(n):
        p = get_parent(i, parent)
        if p not in color:
            stack = [p]
            color[p] = 0
            while stack:
                node = stack.pop()
                for neighbor in adj[node]:
                    if neighbor not in color:
                        color[neighbor] = 1 - color[node]
                        stack.append(neighbor)
                    elif color[neighbor] == color[node]:
                        return False
    return True

def main_solution(n, q, elements):
    elements = [(a-1, b-1, w) for a, b, w in elements]
    
    adj = defaultdict(list)
    parent = list(range(n))
    
    for a, b, w in elements:
        if w == 0:
            set_parent(a, b, parent)
        else:
            if a == b:
                return "no"
            else:
                p1 = get_parent(a, parent)
                p2 = get_parent(b, parent)
                adj[p1].append(p2)
                adj[p2].append(p1)
    
    if is_bipartite(parent, adj, n):
        return "yes"
    else:
        return "no"

# Given input
n = 10
q = 5
elements = [(4, 5, 0), (3, 6, 1), (8, 1, 0), (1, 4, 1), (1, 3, 1)]

# Execute the solution
result = main_solution(n, q, elements)
print(result)